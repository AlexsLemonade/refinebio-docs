# refinebio-docs

## Development workflow

All branches adding or updating documentation should branch from `development`, the default branch for this repository.
Pull requests adding or updating documentation should target `development`.
The development version of the docs can be viewed at <https://docs.refine.bio/en/development/>.

To deploy to latest (<https://docs.refine.bio/en/latest/>), file a pull request to merge `development` into `main`.

### Building on pull requests

Read the Docs builds on pull request events for this repository, which reports the build status and lets you preview changes.

### Spell check

The spell check GitHub Action will fail upon pull requests to `development` if any spelling errors are detected in Markdown files in `docs/`.
The custom dictionary for this repository is located at `config/.custom-dictionary.txt`.

## Local development

### Environment

To set up the environment (assuming `virtualenv` is installed), use the following commands:

```sh
virtualenv env
source env/bin/activate
pip install -r requirements.txt
```

### Building locally

Build locally with the following:

```sh
cd docs/
./autobuild.sh
```

