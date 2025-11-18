"""Utilities for Doxygen documentation generation"""


def doxygenfile_section(
    filename: str,
    project: str,
    project_url: str,
    branch: str = "main",
    dir: str = "src",
) -> str:
    url = f"{project_url}blob/{branch}/{dir}/{filename}"
    ref = filename.replace("_", "-").replace(".", "-")

    # Generate markdown
    return f"""
({ref})=
## `{filename}` <a class="source-link" href="{url}">[source]</a>

```{{doxygenfile}} {filename}
:project: {project}
```
""".strip()


def doxygenfile_sections_toc(
    filenames: list,
    project: str,
    project_url: str,
    branch: str = "main",
    dir: str = "src",
) -> str:
    toc = []
    sec = []
    for filename in filenames:
        ref = filename.replace("_", "-").replace(".", "-")
        toc.append(f"[`{filename}`]({ref})")
        sec.append(
            doxygenfile_section(filename, project, project_url, branch=branch, dir=dir)
        )

    # Create a definition list (": " in markdown) with a newline (" \" in markdown) for each file
    toc = ": " + " \\\n  ".join(toc)

    # Concatenate sections with "---" between each section
    sec = "\n---\n\n".join(sec)

    # Return contents and sections
    return f"""
**Contents**
{toc}

---

{sec}
"""
