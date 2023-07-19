

Looks like pyquante1 has some updating before it can build extension modules under modern Python3.

See:
- https://docs.python.org/3/extending/extending.html
- https://stackoverflow.com/questions/28305731/compiler-cant-find-py-initmodule-is-it-deprecated-and-if-so-what-should-i

I think I'll just generate new bindings using cython rather than try to resurrect the old ones.

Here's how I'm planning to do this:

1. Move cints.c and cints.h to cutil.c and cutil.h.
2. Move wrapper code from cutil.c into cwrap.c
3. The cints_methods list at the end of what is now cwrap.c contains a list of everything we have to wrap in cython. These typically involve creating a cints.pyd file that has a lot of include statements.
4. Would be good to create a test function first. First wrap fact2 (might want to rewrite this along the way).
