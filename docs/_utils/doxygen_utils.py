"""Utilities for Doxygen documentation generation"""


def doxygenfile_section(filename, project, project_url, branch="main", dir="src"):
    url = f"{project_url}blob/{branch}/{dir}/{filename}"
    return f"""
## [`{filename}`]({url})

```{{doxygenfile}} {filename}
:project: {project}
```
""".strip()
