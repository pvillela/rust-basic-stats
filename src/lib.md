A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.

By default, use of this library as a dependency includes all modules. However, each module other than [`core`] has an associated cargo feature that enables the module. To include only selected modules, specify `default-features = false` in the dependency declaration (or `--no-default-features` on the command line) and specify the desired features in the dependency declaration (or command line).
