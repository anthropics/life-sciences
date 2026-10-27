#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Runs one long command in the background and reports on it, for a helper whose bash tool stops any single command
after about 5 minutes (295 seconds). Standard library only. Every file it writes is in the current folder.

    python3 bg.py start NAME -- COMMAND ...   starts COMMAND detached; its output goes to NAME.log, its exit code to NAME.exit
    python3 bg.py wait NAME [SECONDS]         waits up to SECONDS (default 90, at most 120) and prints the job's state
    python3 bg.py status NAME                 prints the job's state at once
    python3 bg.py stop NAME                   kills the job (for example to start it again with changed code)
    python3 bg.py watch NAME [NAME ...]       waits up to 120 seconds (--seconds N first to change it, at most 120)
                                              until none of the jobs is running, and prints all their states

NAME is letters, digits, _ . and -. COMMAND runs without a shell: for a redirect or a pipe, use bash -c '...'.
Each call but watch prints one JSON line: {"ok", "problem", "job", "state", "exit_code", "starts", "log_tail"}. watch,
which the workflow's waiting steps use, prints {"ok", "problem", "state", "jobs"}: state is "running" while any job
runs, else "stopped" if any job stopped, else "done"; jobs gives each job's job, state, exit_code and starts. state is
"running", "done" (the command ended: exit_code 0 is success, another code a failure, a negative one a kill by a signal,
often for lack of memory), or "stopped" (no exit code and no process: the job was killed, for example by a restart of
the machine; start runs it again, and starts counts the starts). A NAME that was never started gives ok false.
"""

import fcntl
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path


def files(name):
    return {k: Path(name + "." + k) for k in ("log", "exit", "pid", "starts", "lock")}


def alive(name):
    """True while the job's runner process lives. The runner holds a lock on NAME.lock from before it starts until it
    ends, however it ends, so the job is alive exactly when that lock cannot be taken."""
    try:
        with open(files(name)["lock"], "rb") as h:
            fcntl.flock(h, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return True
    except OSError:  # no lock file, or a file system without locks, on which start refuses to start a job
        pass
    return False


def state(name):
    f = files(name)
    starts = int(f["starts"].read_text()) if f["starts"].is_file() else 0
    tail = ""
    if f["log"].is_file():
        with open(f["log"], "rb") as h:
            h.seek(max(0, f["log"].stat().st_size - 600))
            tail = h.read().decode("utf-8", "replace")
    if f["exit"].is_file():
        try:
            code = int(f["exit"].read_text().strip())
        except ValueError:
            code = None
        return {"job": name, "state": "done", "exit_code": code, "starts": starts, "log_tail": tail}
    if not starts:
        return {"job": name, "state": "not started", "exit_code": None, "starts": 0, "log_tail": ""}
    running = alive(name)
    if not running and f["pid"].is_file() and time.time() - f["pid"].stat().st_mtime < 10:
        time.sleep(0.5)  # a job started a moment ago can read as gone for an instant while its process starts: look again
        running = alive(name)
    if running:
        return {"job": name, "state": "running", "exit_code": None, "starts": starts, "log_tail": tail}
    if f["exit"].is_file():  # the job ended between the two checks above
        return state(name)
    return {"job": name, "state": "stopped", "exit_code": None, "starts": starts, "log_tail": tail}


def watch(a, usage):
    """Waits until none of the named jobs is running (at most 120 s) and prints their states in one JSON line."""
    limit = 120.0
    if len(a) >= 2 and a[0] == "--seconds" and re.fullmatch(r"[0-9]+(\.[0-9]+)?", a[1]):
        limit, a = min(120.0, float(a[1])), a[2:]
    if not a or not all(re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", n) for n in a):
        print(json.dumps({"ok": False, "problem": usage}))
        sys.exit(2)
    unknown = [n for n in a if state(n)["state"] == "not started"]
    if unknown:
        print(json.dumps({"ok": False, "problem": "no job named %s in this folder" % ", ".join(unknown), "state": "", "jobs": []}))
        return
    end = time.time() + limit
    while time.time() < end and any(state(n)["state"] == "running" for n in a):
        time.sleep(min(2.0, max(0.1, end - time.time())))
    jobs = [{k: s[k] for k in ("job", "state", "exit_code", "starts")} for s in map(state, a)]
    states = {j["state"] for j in jobs}
    overall = "running" if "running" in states else "stopped" if "stopped" in states else "done"
    print(json.dumps({"ok": True, "problem": "", "state": overall, "jobs": jobs}))


def main():
    a = sys.argv[1:]
    if len(a) >= 2 and a[0] == "_run":
        # the detached runner: runs the command, then writes its exit code last (so NAME.exit means it ended)
        name, cmd = a[1], a[3:]
        with open(name + ".log", "ab") as log:
            try:  # unbuffered, so a Python job's progress and errors reach NAME.log as they happen
                code = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                       env=dict(os.environ, PYTHONUNBUFFERED="1"))
            except OSError as e:  # the command could not be started at all (for example a missing python3 in front)
                log.write(("bg.py: could not start the command: %s\n" % e).encode())
                code = 127
        Path(name + ".exit.tmp").write_text(str(code))
        os.replace(name + ".exit.tmp", name + ".exit")
        return
    usage = "usage: bg.py start NAME -- COMMAND ... | wait NAME [SECONDS] | status NAME | stop NAME | watch [--seconds N] NAME [NAME ...] (NAME: letters, digits, _ . -)"
    if a and a[0] == "watch":
        watch(a[1:], usage)
        return
    ok_args = (len(a) >= 2 and a[0] in ("start", "wait", "status", "stop") and re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", a[1])
               and (a[0] != "start" or (len(a) >= 4 and a[2] == "--")) and (a[0] != "wait" or len(a) < 3 or re.fullmatch(r"[0-9]+(\.[0-9]+)?", a[2])))
    if not ok_args:
        print(json.dumps({"ok": False, "problem": usage}))
        sys.exit(2)
    op, name = a[0], a[1]
    f = files(name)
    if op != "start" and state(name)["state"] == "not started":
        print(json.dumps(dict(state(name), ok=False, problem="no job named %s in this folder: start it first" % name)))
        return
    if op == "stop":
        if state(name)["state"] == "running":
            with open(f["log"], "ab") as log:
                log.write(b"\n--- stopped by bg.py stop\n")
            try:
                os.killpg(int(f["pid"].read_text().strip()), signal.SIGKILL)  # the runner leads its own process group
            except (OSError, ValueError):
                pass
            Path(name + ".exit").write_text(str(-signal.SIGKILL))
        print(json.dumps(dict(state(name), ok=True, problem="")))
        return
    if op == "start":
        if state(name)["state"] == "running":
            print(json.dumps(dict(state(name), ok=False, problem="already running: wait for it, or stop it first if you changed its code")))
            return
        with open(f["lock"], "wb") as lock:  # taken here and inherited by the runner, which keeps it when this process ends
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print(json.dumps(dict(state(name), ok=False, problem="%s is being started or checked by another command right now: run bg.py status %s, and start it again only if it is not running" % (name, name))))
                return
            except OSError as e:
                print(json.dumps(dict(state(name), ok=False, problem="could not lock %s.lock (%s): run in a folder on a local disk" % (name, e))))
                return
            for k in ("exit", "pid"):
                f[k].unlink(missing_ok=True)
            starts = (int(f["starts"].read_text()) if f["starts"].is_file() else 0) + 1
            f["starts"].write_text(str(starts))
            with open(f["log"], "ab") as log:
                log.write(("\n--- start %d: %s\n" % (starts, " ".join(a[3:]))).encode())
            p = subprocess.Popen([sys.executable, os.path.abspath(__file__), "_run", name] + a[2:], start_new_session=True, pass_fds=(lock.fileno(),),
                                 stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        f["pid"].write_text(str(p.pid))
        print(json.dumps(dict(state(name), ok=True, problem="")))
        return
    if op == "wait":
        limit = min(120.0, float(a[2])) if len(a) > 2 else 90.0
        end = time.time() + limit
        while time.time() < end and state(name)["state"] == "running":
            time.sleep(2)
    print(json.dumps(dict(state(name), ok=True, problem="")))


if __name__ == "__main__":
    main()
