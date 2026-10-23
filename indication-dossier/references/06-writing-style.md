# Writing Style

Write in the style of an industry-trained analyst presenting the facts clearly and completely.

# Citation Standards

## Inline Citation Format

**Use inline markdown links for all factual claims:**

```markdown
[descriptive claim text](source_url)
```

**Examples:**
- "Efgartigimod achieved [68% MG-ADL responder rate vs 30% placebo](https://pubmed.ncbi.nlm.nih.gov/34146511/)"
- "Phase 3 [ADAPT trial (NCT03669588)](https://clinicaltrials.gov/study/NCT03669588)"

## When to Use Inline Links vs Reference Lists

**Prefer inline markdown links** for:
- Quantitative data that readers may want to verify immediately
- Clinical trial results and endpoints
- Regulatory approvals and dates
- Competitive claims about stage, timing, efficacy

**Use numbered reference lists [1], [2]** only when:
- The same source is cited 5+ times (define once, reference by number)
- Citing multiple sources for a single claim
- The URL is excessively long (>150 characters) and would break reading flow

## URL Best Practices

**Link specificity:**
- ClinicalTrials.gov: Link directly to study page with NCT number
  - Good: `https://clinicaltrials.gov/study/NCT03669588`
  - Bad: `https://clinicaltrials.gov` (homepage)

- PubMed: Use PubMed ID link format
  - Good: `https://pubmed.ncbi.nlm.nih.gov/34146511/`
  - Avoid: DOI redirects when PubMed ID available

- Company materials: Link to specific press release or presentation
  - Good: `https://www.company.example/news/company-announces-phase-3-results`
  - Bad: `https://www.company.example` (requires navigation)

- Patents: Use Google Patents for accessibility
  - Good: `https://patents.google.com/patent/US########X#`
  - Alternative: USPTO link if needed for official text

## Citation Discipline: When to Cite

**ALWAYS cite (with inline link):**
- Quantitative data (patient numbers, efficacy metrics)
- Clinical endpoints and results (response rates, mean changes, p-values)
- Development status (who's developing what, at what stage)
- Dates and timelines (approval years, trial completion estimates)
- Patent numbers and expiry dates
- Safety data (adverse event rates, discontinuation rates, black box warnings)

**DON'T cite:**
- General medical knowledge (e.g., "Phase 3 trials are larger than Phase 2")
- Your own analytical interpretations
- Mathematical calculations you perform from cited data (cite the inputs, not the calculation)

## Citation Quality Hierarchy

Follow the source priority in `00-research-standards.md`, so research and writing resolve conflicts the same way.

## Handling Conflicting Sources

When sources disagree, cite both and explain:

```markdown
The press release reports [200 patients](source1), while ClinicalTrials.gov shows [180 enrolled](source2).
The registry figure of 180 is used.
```

## Clinical Trial Citations

**Format for trial references:**
- Include NCT number as a clickable link
- Add trial name/acronym if available
- Example: "Phase 3 [ADAPT trial (NCT03669588)](https://clinicaltrials.gov/study/NCT03669588)"

**When citing trial results:**
- Link to publication if published
- Link to ClinicalTrials.gov if results posted there
- Link to company press release if that's the only source
- Example: "[ADAPT trial](https://pubmed.ncbi.nlm.nih.gov/34146511/) demonstrated 68% MG-ADL responder rate"

## Quality Checks Before Finalizing

Before submitting your analysis, verify:

- [ ] Every quantitative claim has an inline citation
- [ ] Every therapy in development has stage + trial ID cited
- [ ] All approval dates and timelines are cited
- [ ] Clinical endpoints and efficacy data link to trials or publications
- [ ] URLs are specific (not homepages) and functional
- [ ] Conflicting sources are acknowledged and resolved

## Common Citation Mistakes to Avoid

**DON'T:**
- "Studies show..." (which studies? cite them)
- "Drug X is in Phase 3" (link to trial registry)
- "The drug was approved in 2023" (link to FDA announcement or label)
- Citing Wikipedia, discussion boards, or unreliable sources
- Backfilling citations after writing (cite as you write)

**DO:**
- "[Drug X's Phase 3 trial (NCT#####)](url) enrolled 200 patients"
- "[FDA approved the drug in June 2023](url) for adult patients"
- Add citations immediately when you extract information from sources
