# TODOs

## Getting started with Julia

> :info: Make sure you have the `dot` program installed on your computer to be able to visualize graphs.

1. Download and install the latest version (1.9.1): [https://julialang.org/downloads/]()
2. Using Jupyter notebook (optional but recommended)
```bash
$ julia
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.1 (2023-06-07)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
(@v1.9) pkg>  # <- enter the package mode by hitting the "]"
(@v1.9) pkg> add Pluto

# After the installation
julia> using Pluto
julia> Pluto.run() # will open a new tab in your browser.
```
3. To develop a package and debug in the notebook:
  1. Clone the repository
  2. in the first cell of the notebook:
  ```julia
  begin
      using Pkg

      # Path to the clone repository.
      # you need to save the notebook to know the relative path
      Pkg.develop(path="...")

      # Take into account the changes you've made
      using Revise

      using FiniteStateAutomata # or Semirings or SparseSemimodules
  end
  ```

## Semirings [Semirings.jl](https://gitlab.lisn.upsaclay.fr/fast/semirings.jl)
- [ ] Semirings
  - [x] Boolean 🌶
  - [x] Real 🌶
  - [x] LogSemiring 🌶
  - [x] Tropical (max-plus) 🌶
  - [ ] Tropical (min-plus) 🌶
  - [ ] Expectation 🌶
  - [ ] Lexicographic 🌶
  - [ ] MinMax 🌶
  - [ ] Power 🌶
  - [x] Product 🌶
  - [ ] SignedLog 🌶
  - [ ] SparsePower 🌶
  - [ ] String 🌶
- [ ] Automatic differentiation (two possible choices)
  - [ ] [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) compute the
    gradient from the LLVM IR 🌶🌶🌶
  - [ ] [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl)
    with [Zygote.jl](https://github.com/FluxML/Zygote.jl) source-to-source
    differentiation 🌶🌶🌶

## Linear Algebra [SparseSemimodules.jl](https://gitlab.lisn.upsaclay.fr/fast/sparsesemimodules.jl)
- [x] Basic structures and algorithms
  - [x] sparse vector
  - [x] sparse matrix CSR
  - [x] dot product of sparse vectors
  - [x] sparse vector transpose - sparse matrix CSR multiplication
  - [ ] kronecker Product (for the composition algorithm)
    - [ ] Kronecker product of sparse matrix CSR 🌶🌶
    - [ ] Lazy version of the Kronecker product 🌶🌶
- [ ] Parallel computations
  - [ ] Sum of sparse matrices 🌶🌶
  - [ ] Multi-threaded sparse vector tr. - sparse matrix multiplication
    🌶🌶
  - [ ] GPU sparse vector tr. - sparse matrix multiplication 🌶🌶🌶
- [ ] Matrix star operations 🌶
  - [x] example [here](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl/-/blob/jsalt2023-workshop/src/transitionmatrix.jl)
    to be relocated in SparseSemimodules.jl 🌶
  - [ ] Option to prune the paths 🌶🌶

## FST [FiniteStateAutomata.jl](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl)
- [ ] ⚠ **FST storage format** ⚠ (high priority)
  - [x] List of sparse matrices CSR stored as a block diagonal matrix
  - [ ] Other alternatives ?
  - [ ] efficient arcs / state iterators
- [ ] Visualization
  - [x] graphviz
  - [ ] show summary
- [ ] Automata's weight (shortest distance) 🌶
  - [x] basic algorithm
  - [ ] batch version
- [ ] Rational operations
  - [ ] union 🌶
  - [ ] concatenation 🌶🌶
  - [ ] closure 🌶🌶
- [ ] Composition / Intersection
  - [ ] ϵ-free 🌶🌶
  - [ ] with ϵ arcs 🌶🌶🌶
  - [ ] lazy implementation 🌶 (relies on SparseSemimodules)
- [ ] Weight pushing 🌶🌶
- [ ] Other FST Operations
  - [ ] determinize 🌶🌶🌶
  - [ ] minimize 🌶🌶🌶
  - [ ] ϵ-arc removal 🌶🌶
  - [ ] reverse 🌶
  - [ ] invert 🌶
  - [ ] project 🌶
  - [ ] connect 🌶
    - [ ] requires to implement: find accessible / coaccessible states 🌶
  - ...
- [ ] Python wrapper 🌶🌶
- [ ] OpenFST bindings
  - [x] OpenFST text format 🌶
  - [ ] "binary binding" (Michael ?)

## Practical uses
- [ ] Decoder (edge-cloud / long-form pipelines ?)
- [ ] CTC / MMI Loss
  - [ ] Create a language model from a IARPA n-gram file
  - [ ] requires local and global normalization (i.e. weight pushing)
  - [ ] Intersection algorithm for GNAT ??
- [ ] SoundsLikeName
  - [ ] pronunciation graph for OOV

