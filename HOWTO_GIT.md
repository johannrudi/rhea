# HowTo - Git Guidelines

## Commits

### Commit messages

Recommended guidelines:
- Start sentences/messages with a capital letter.
- Start sentences/messages with a verb.
- End the commit message with a "."
- If commiting files from `src` directory, then add the filenames without
  extension, for example, `git commit -m "rhea_viscosity: Did this, because of that."`
- If commiting files from an example (e.g., the `basic` example), then start
  the commit message with "Example basic: Did this, because of that."
- A commit message should talk about **what** changed, and **why**. Not
  **how**, because the diff shows how and it does not need to be repeated.

### Commit workflow

Break up code changes into (small) incremental commits. Use `git add ...` and
`git commit -m "..."` workflow to create commits of specific files.

#### Example:
Assume we changed the the following files:
- `rhea_viscosity.c` and `rhea_viscosity.h` in the `src` directory (e.g.,
  fixing a bug in the Arrhenius law)
- `basic.c` in the `example/basic` directory (e.g., adding output in HDF5 format)

Then the recommended way of commiting these changes is:
```bash
cd src
git add rhea_viscosity.c rhea_viscosity.h
git commit -m "rhea_viscosity: Fix bug in Arrhenius law, because it caused NAN values."

cd ../example/basic
git add basic.c
git commit -m "Example basic: Add output in HDF5 format."
```

## Branches

Recommended:
- Use branches that are created from issues (on Bitbucket's web interface).
- When finished working on one branch, start a pull request (on Bitbucket's web interface).

Alternative:
- Create your custom developemnt branch, starting the branch name with
  `your name + slash + develop`:
  ```bash
  git branch mike/develop
  git checkout mike/develop
  ```

## Fetching/Pulling

Recommended: Use `git fetch origin` then `git merge origin/some-branch` instead
of `git pull`!

#### Example:
```bash
git fetch origin
# ... prints which remote branches have been updated/created since you last fetched

git branch -vv
# ... prints how local branches are ahead/behind remote branches

git merge origin/some-branch
# ... merges the remote branch "some-branch" into the currently checked out local branch
```

## Pushing

Recommended: Use specific push commands that specify the remote name and the branch name!

#### Example:
Assume you want to push your custom development branch, which is called
`mike/develop` to the remote called `origin`; then run:
```bash
git push origin mike/develop
```

