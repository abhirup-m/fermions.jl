# Fermions.jl 

*A Toolkit for Working with Models of Interacting Electrons*
 
![](./examples/images/sshEdge.svg)

**Fermions.jl** is a toolkit for designing and analysing second-quantised many-particle Hamiltonians of electrons, potentially interacting with each other. The main point in designing this library is to abstract away the detailed task of writing matrices for many-body Hamiltonians and operators (for correlations functions) with large Hilbert spaces; all operators (including Hamiltonians) can be specified using predefined symbols, and the library then provides functions for diagonalising such Hamiltonians and computing observables within the states. (In case you are not accustomed to using second-quantised operators, check [this brief explanation](#a-brief-explanation-of-second-quantised-operators-for-the-uninitiated).)

## Neat features

- Construct many-body operators and states using simple datastructures such as dictionaries and vectors. No need to bother with complicated and abstract classes and their objects.

- High-level of freedom in constructing fermionic Hamiltonians. All Hamiltonians that can be represented as a tensor product of 2-dimensional fermionic Fock-space operators can be modelled using fermions.jl.

- Uses optimised algorithms that make use of symmetries of the problem. Models with total occupancy conservation and/or total magnetisation conservation can be automatically block-diagonalised.

- Provides a wide range of inbuilt functions for calculating various quantities of physical interest, including spectral functions, static correlations and measures of entanglement. The ability to construct any general correlation function by using fermionic operators further extends the range of possibilities.

This library was borne out of a need to numerically construct and solve fermionic Hamiltonians in the course of my doctoral research. While there are similar julia libraries such as [Marco-Di-Tullio/Fermionic.jl](https://github.com/Marco-Di-Tullio/Fermionic.jl) and [qojulia/QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl), **fermions.jl is much more intuitive** since it works directly on predefined basis states and allows defining arbitrary fermionic operators and quantum mechanical states. **There is no need to interact with complicated and abstract classes and objects** in order to use this library; everything is defined purely in terms of simple datastructures such as dictionaries, vectors and tuples. This makes the entire process transparent and intuitive.

## Will this be useful for me?

You might find this library useful if you spend a lot of time studying Hamiltonian models of fermionic or spin-1/2 systems, particularly ones that cannot be solved analytically, or use a similar library in another language (QuTip in python, for example), but want to migrate to Julia. You will not find this useful if you mostly work with bosonic systems and open quantum systems, or work in the thermodynamic limit (using methods like quantum Monte Carlo, numerical RG).

## Documentation Outline
```@contents
Pages = ["quickstart.md", "base.md", "eigen.md", "correlations.md", "theory.md"]
Depth = 1
```

Feedback and contributions are always welcome!
