# Fermions.jl - Create and Analyse Models of Interacting Electrons

## What is this?
**Fermions.jl** is a Julia library for designing and analysing second-quantised many-particle Hamiltonians of electrons, potentially interacting with each other. The main point in designing this library is to abstract away the detailed task of writing matrices for many-body Hamiltonians and operators (for correlations functions) with exponentially large sizes; all operators (including Hamiltonians) can be specified using symbols, and the library then provides functions for diagonalising such Hamiltonians and computing observables within the states.
<details>
  <summary><b>What are second-quantised operators?</b></summary>
  Second quantised operators are a convenient way of writing many-body Hamiltonians. As an example, we wish to design the Hamiltonian of a single electron hopping across a 1D chain of lattice sites labelled as $i=1$, $i=2$, $i=3$ and so on. The system is such that the electron can only hop on to the nearest-neighbour sites starting from any given site: $i\to i+1$, $i\to i-1$. In order to hop from, say, $i=1$ to $i=2$, the electron must first be _annihilated_ at the site $i=1$ and then _created_ at $i=2$. The creation process is represented by the operator $c^\dagger_1$ (1 representing the site index of the location of the operation), while the annihilation process is represented by $c_2$, the total process being the product of both process: $c^\dagger_1 c_2$. Of course, the opposite process can also happen - the electron can also hop from $i=2$ to $i=1$, so that the total second-quantised Hamiltonian for the dynamics involving sites 1 and 2 is 
$$
c^\dagger_1 c_2 + c^\dagger_2 c_1~.
$$
Now, the indices 1 and 2 can be represented by any consecutive indices $i$ and $i+1$, leading to the so-called tight-binding Hamiltonian in 1-dimensions:
$$
H_\mathrm{TB} = \sum_{i} \left(c^\dagger_i c_{i+1} + c^\dagger_{i+1} c_i\right)
$$
</details>
