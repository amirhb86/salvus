# Build salvus' Documentation

`salvus'` documentation is build by employing a combination of `Doxygen` (for
auto-generated C++ documentation) and `sphinx` (to render stuff). The `breathe`
project ties both together.

You must have a recent version of `Doxygen` installed, as well as some Python
version. Install the Python dependencies and build the documentation with:

```bash
$ pip install sphinx sphinx_rtd_theme breathe
$ cd /path/to/salvus
$ cd doc
$ make html
```

Now just open `_build/html/index.html`.

The `make html` command runs `Doxygen` and will auto-generate a `Doxygen`-like
documentation which is then rendered with `sphinx`.
