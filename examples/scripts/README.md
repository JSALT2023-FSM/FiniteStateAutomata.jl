
# Examples

A few examples of how to use the toolkit to create and manipulate FSTs.

## Visualization

For **small FSTs** it is possible to display the graphical representation
directly in the terminal. Supported terminals are:
- [iTerm2](https://iterm2.com) (MacOS)
- [xterm](https://en.wikipedia.org/wiki/Xterm) (linux)

!!! info
   On linux you need to start `xterm` with the following command:
    ```
    xterm -fa 'Monospace' -fs 14 -xrm "XTerm*decTerminalID: vt340" -xrm "XTerm*numColorRegisters: 256"
    ```

## Running the examples

First, instantiate the virtual environment:
```
$ julia --project=./ -e 'using Pkg; Pkg.instantiate()`
```
then run the example you want, for instance:
```
$ julia --project=./ 1_createfst.jl
```

