# TODOs

## Semirings [Semirings.jl](https://gitlab.lisn.upsaclay.fr/fast/semirings.jl)
- Semiring to implement
  - [x] Boolean 🌶
  - [x] Real 🌶
  - [x] LogSemiring 🌶
  - [x] Tropical (max-plus) 🌶
  - [ ] Tropical (min-plus) 🌶
  - [ ] Expectation 🌶
  - [ ] Lexicographic 🌶
  - [ ] MinMax 🌶
  - [ ] Power 🌶
  - [ ] Product 🌶
  - [ ] SignedLog 🌶
  - [ ] SparsePower 🌶
  - [ ] String 🌶
- Automatic differentiation (two possible choices)
  - [ ] [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) compute the
    gradient from the LLVM IR 🌶🌶🌶
  - [ ] [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl)
    with [Zygote.jl](https://github.com/FluxML/Zygote.jl) source-to-source
    differentiation 🌶🌶🌶

## Linear Algebra [SparseSemimodules.jl](https://gitlab.lisn.upsaclay.fr/fast/sparsesemimodules.jl)
- [ ] Parallel computations
  - [x] Sum of sparse matrices 🌶🌶
  - [ ] Multi-threaded sparse vector - sparse matrix multiplication
    🌶🌶
  - [ ] GPU sparse vector - sparse matrix multiplication 🌶🌶🌶
- [ ] Kronecker Product (for the composition algorithm)
  - [ ] Kronecker product of sparse matrix CSR 🌶🌶
  - [ ] Lazy version of the Kronecker product 🌶🌶
- [ ] Matrix star operations 🌶
  - [ ] Example [here](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl/-/blob/jsalt2023-workshop/src/transitionmatrix.jl)
    to be relocated in SparseSemimodules.jl 🌶
  - [ ] Option to prune the paths 🌶

## FST [FiniteStateAutomata.jl](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl)
- [ ] **FST storage format**
  - [x] List of sparse matrices CSR stored as a block diagonal matrix
  - [ ] Other alternatives ?
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
- [ ] OpenFST binding 🌶🌶 (Michael Riley ?)

