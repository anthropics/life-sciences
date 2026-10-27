"""Tests for check_repo.py. Run with: python3 -m unittest discover -s .github/scripts"""

import json
import tempfile
import unittest
from pathlib import Path

import check_repo

SKILL = "---\nname: good-skill\ndescription: Does a thing.\n---\n\n# Good skill\n"


class CheckRepoTest(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.write("good-skill/SKILL.md", SKILL)
        self.write("good-skill/LICENSE.txt", "Apache 2.0")
        self.write("README.md", "# Title\n\n[Skill](good-skill/SKILL.md)\n\n```bash\nls\n```\n")
        marketplace = {"plugins": [{"name": "good-skill", "skills": ["./good-skill"]}]}
        self.write(".claude-plugin/marketplace.json", json.dumps(marketplace))

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def problems(self):
        return (
            check_repo.check_marketplace(self.root)
            + check_repo.check_markdown(self.root)
            + check_repo.check_data_files(self.root)
        )

    def assert_finds(self, text):
        self.assertTrue(any(text in problem for problem in self.problems()), self.problems())

    def test_valid_repository_passes(self):
        self.assertEqual(self.problems(), [])

    def test_finds_unclosed_code_block(self):
        self.write("README.md", "# Title\n\n```bash\nls\n")
        self.assert_finds("README.md:3: code block is never closed")

    def test_longer_fence_can_hold_a_shorter_one(self):
        self.write("README.md", "````\n```\ninside\n```\n````\n")
        self.assertEqual(self.problems(), [])

    def test_finds_broken_link(self):
        self.write("README.md", "[Gone](missing/file.md)\n")
        self.assert_finds("broken link to missing/file.md")

    def test_ignores_links_in_code_web_links_and_anchors(self):
        self.write("README.md", "```\n[x](nope.md)\n```\n`[y](nope.md)` [web](https://a.b) [top](#title)\n")
        self.assertEqual(self.problems(), [])

    def test_ignores_format_example_links(self):
        self.write("README.md", "Cite as [title](url).\n")
        self.assertEqual(self.problems(), [])

    def test_finds_invalid_json_and_yaml(self):
        self.write("good-skill/data.json", "{")
        self.write("good-skill/data.yaml", "a: [")
        self.assert_finds("invalid JSON")
        self.assert_finds("invalid YAML")

    def test_reads_yaml_with_several_documents(self):
        self.write("good-skill/two.yaml", "a: 1\n---\nb: 2\n")
        self.assertEqual(self.problems(), [])

    def test_finds_listed_skill_without_skill_file(self):
        (self.root / "good-skill/SKILL.md").unlink()
        self.assert_finds("has no SKILL.md")

    def test_finds_skill_missing_from_marketplace(self):
        self.write("stray/SKILL.md", SKILL)
        self.assert_finds("stray: has a SKILL.md but is not in the marketplace")

    def test_finds_bad_frontmatter(self):
        self.write("good-skill/SKILL.md", "---\nname: Bad Name\n---\n")
        self.assert_finds("name must be")
        self.assert_finds("description is missing")

    def test_finds_missing_license(self):
        (self.root / "good-skill/LICENSE.txt").unlink()
        self.assert_finds("LICENSE.txt is missing")


if __name__ == "__main__":
    unittest.main()
