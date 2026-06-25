# Releasing

## GitHub

Push changes to `master` and let the GitHub Actions test workflow pass.

## PyPI

This repository includes a publish workflow at `.github/workflows/publish.yml`.

To finish the setup, one of the following must be configured in GitHub/PyPI:

1. Trusted publishing:
   - Create the project on PyPI if it does not exist.
   - In PyPI, add this GitHub repository as a trusted publisher.
   - In GitHub, keep the `pypi` environment available for the publish workflow.
   - Recommended owner/repository: `davamaro/stabgraph`.

2. Manual publish from a local machine:
   - Install `build` and `twine`.
   - Run `python -m build`.
   - Run `python -m twine check dist/*`.
   - Run `python -m twine upload dist/*`.

## Suggested release flow

1. Update `CHANGELOG.md`.
2. Bump the version in `setup.py`.
3. Build locally with `python -m build --no-isolation` if your current machine
   has restricted network access.
4. Validate with `python -m twine check dist/*`.
5. Create a GitHub release or run the `publish` workflow manually.
