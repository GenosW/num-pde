# Numerical methods for PDEs

Lecturer: <https://www.asc.tuwien.ac.at/~schoeberl/wiki/index.php/Philip_lederer>

Homepage: <https://www.asc.tuwien.ac.at/~schoeberl/wiki/index.php/Numpde20>

NGSolve: <https://ngsolve.org/>

Lectures:   Mo., 14:00 - 16:00
            Tu., 14:00 - 15:00

Exercises:  Tu., 15:00 - 17:00

---

## NGSolve/netgen

NGSolve-6.2.1910 --> python 3.7

newer --> python 3.8 (NO!)

For MacOS X using python managed by homebrew:

```bash
sudo ln -s /usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7m.dylib /Library/Frameworks/Python.framework/Versions/3.7/Python
```

The python DYLIB needs to be symlinked to the framework folder/path expected by NGSolve.
Adjust accordingly for different version of python.

Also, set some env variables

```bash
export PYTHONPATH=$PYTHONPATH:/Applications/Netgen.app/Contents/Resources/lib/python3.7/site-packages:.
export NETGENDIR=/Applications/Netgen.app/Contents/MacOS
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$NETGENDIR
export DYLD_FRAMEWORK_PATH=$DYLD_FRAMEWORK_PATH:$NETGENDIR/../Frameworks
export PATH=$NETGENDIR:$PATH
```

## Windows

Cmder alias: 

```bash
alias A-pde=D:\Studium\Code\num-pde\.env\Scripts\activate $T cd D:\Studium\Code\num-pde
```

[Similiar issue/solution](https://ngsolve.org/forum/ngspy-forum/163-macos-installation-frameworks-path-missing-at-ngsolve-6-2-1804-dmg)

## 1. Lecture

## 2. Lecture
