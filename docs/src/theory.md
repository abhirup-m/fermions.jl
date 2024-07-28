# Theoretical Background - Fermions, Symmetries and Second Quantisation

This page lays out some of the theoretical details that have been utilised while developing this toolkit. All these ideas are very rudimentary from the perspective of practising condensed matter physicists, so this section is mostly for the sake of completeness.

## A Brief explanation of second-quantised operators
The entire library is built on the idea of using second-quantisation to study many-body Hamiltonians, so a discussion of that is certainly pertinent. Second quantised operators are a convenient way of writing many-body Hamiltonians. This library is only concerned with electrons; these are particles whose states must be antisymmetric with respect to an exchange of the particles. For example, a state of electrons at positions 1 and 2 has to be $|1,2\rangle - |2, 1\rangle$. In order to avoid this clumsy notation, we introduce creation and annihilation operators $c^\dagger_\nu$ and $c_\nu$, which carry the property of antisymmetry within themselves. The two operators, respectively, create and annihilate an electron in the state described by the quantum number $\nu$, which can be the position or momentum, among others; they satisfy $c_1 c_2 = -c_2 c_1$, allowing us to write the above complicated state simply as $c^\dagger_1 c^\dagger_2|0\rangle$.

As an example, we wish to design the Hamiltonian of a single electron hopping across a 1D chain of lattice sites labelled as $i=1$, $i=2$, $i=3$ and so on. The system is such that the electron can only hop on to the nearest-neighbour sites starting from any given site: $i\to i+1$, $i\to i-1$. In order to hop from, say, $i=1$ to $i=2$, the electron must first be _annihilated_ at the site $i=1$ and then _created_ at $i=2$. The creation process is represented by the operator $c^\dagger_1$ (1 representing the site index of the location of the operation), while the annihilation process is represented by $c_2$, the total process being the product of both process: $c^\dagger_1 c_2$. Of course, the opposite process can also happen - the electron can also hop from $i=2$ to $i=1$, so that the total second-quantised Hamiltonian for the dynamics involving sites 1 and 2 is 

$$c^\dagger_1 c_2 + c^\dagger_2 c_1~.$$

Now, the indices 1 and 2 can be represented by any consecutive indices $i$ and $i+1$, leading to the so-called tight-binding Hamiltonian in 1-dimensions:

$$H_\mathrm{TB} = \sum_{i} \left(c^\dagger_i c_{i+1} + c^\dagger_{i+1} c_i\right)$$
