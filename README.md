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
    [Zygote.jl](https://github.com/FluxML/Zygote.jl) source-to-source
    differentiation 🌶🌶🌶

## Linear Algebra [SparseSemimodules.jl](https://gitlab.lisn.upsaclay.fr/fast/sparsesemimodules.jl)
- [ ] Parallel computations
  - [x] Sum of sparse matrices 🌶🌶
  - [ ] Multi-threaded sparse vector - sparse matrix multiplication
    🌶🌶
  - [ ] GPU sparse vector - sparse matrix multiplication 🌶🌶🌶
- [-] Kronecker Product (for the composition algorithm)
  - [-] Kronecker product of sparse matrix CSR 🌶🌶
  - [-] Lazy version of the Kronecker product 🌶🌶
- [-] Matrix star operations 🌶
  - [-] Example [here](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl/-/blob/jsalt2023-workshop/src/transitionmatrix.jl)
    to be relocated in SparseSemimodules.jl
  - [ ] Option to prune the paths

## FST [FiniteStateAutomata.jl](https://gitlab.lisn.upsaclay.fr/fast/finitestateautomata.jl)
- [-] Automata's weight (shortest distance) 🌶
- [ ] Connect
  - [ ] find accessible / coaccessible states
  - [ ]
- [ ] Composition
  - [ ] ϵ-free
  - [ ] with ϵ arcs
  - [ ] lazy algorithm
- [ ] Composition
  - [ ] ϵ-free
  - [ ] with ϵ arcs
  - [ ] lazy algorithm
- [ ] Python wrapper 🌶🌶
  -
