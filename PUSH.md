# Publishing the artifact

The E&D track requires the data and code to be reachable by reviewers **at submission time**. The default review
mode is double-blind, but the call for papers allows declaring single-blind when the artifact cannot be anonymised.

**Decision taken for this submission: single-blind, hosted at
<https://github.com/tienenCUHK/developability-bench>.** Declare single-blind on the submission form, and keep the
repository public for the duration of review.

Everything needed is prepared: `release/` is a git repository with two commits, author
`Anonymous Authors <anonymous@example.org>`, no usernames, no server paths, no acknowledgements, and a checksum
file. A portable copy also exists as `devmpnn-bench-anon.bundle`.

If the decision is ever reversed, the three standard routes are, in order of preference:

## Route A — anonymous GitHub account (most common)

1. Create a fresh GitHub account with a neutral name (no institutional e-mail, no personal handle).
2. Create a **public** repository, e.g. `developability-bench`.
3. From this directory:

   ```bash
   cd release
   git remote add origin https://github.com/<neutral-account>/developability-bench.git
   git branch -M main
   git push -u origin main
   ```

4. Sanity-check the rendered README in a private browser window; nothing in it should be linkable to you.
5. Put the URL in the submission form. Update or delete the account after the review period if you prefer.

## Route B — anonymous mirror of an existing repository

Services such as `anonymous.4open.science` create an anonymised mirror of a GitHub repository for the duration of
review. They need a source repository first, so this is Route A plus one extra step: paste the repository URL into
the mirror service and submit the mirror link instead. Use this if you would rather not maintain a second account
for long.

## Route C — private link (last resort)

Some authors share a read-only archive link (e.g. a Zenodo record created under a pseudonym, or a file drop).
A bare Google Drive or personal-institution link is a common desk-reject trigger for this track because it can be
unreachable or identity-revealing, so prefer A or B.

## Data hosting (needed in addition to the code repository)

The CFP asks for datasets on a dedicated host (Hugging Face, Dataverse, Kaggle, OpenML) and requires a Croissant
file with core **and** RAI fields.

- Preferred hosts auto-generate the Croissant core file. If you use one, merge the RAI block from
  `croissant.json` into the generated file rather than replacing it wholesale, then validate.
- Validate with the official checker before submitting; the file in this repository has been checked
  field-by-field against the Croissant RAI 1.0 property table, but the online checker is the authority.

Upload targets, in order of size:

| target | content | size |
|---|---|---|
| dataset | `data/b30_designs.tsv`, `data/p16/`, `data/lam_sweep/`, `data/px_scTM.tsv`, `data/targets30.txt` | ~2.6 MB |
| code | `code/`, `paper/`, `DATASHEET.md`, `BENCHMARK_CARD.md`, `croissant.json` | ~1.3 MB |

Both fit comfortably in a free Hugging Face dataset repository.

## After uploading

1. Put both links in `SUBMISSION.md` and in the OpenReview submission form.
2. Re-run `shasum -a 256 -c CHECKSUMS.sha256` inside `release/` to confirm the uploaded tree matches what the
   paper was built from.
3. If anonymity later turns out to be broken (for example a personal e-mail in a git commit), rewrite the history
   rather than deleting the repository: `git commit --amend --author` on both commits, then force-push.
