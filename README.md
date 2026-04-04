# truss-solver
A header-only C++ project for solving internal forces and nodal displacements in 2D/3D truss systems using Method of Joints (Static Equilibirum) and Finite Element Method. 


## Derivation of Stiffness Matrix $[k_e]$

Displacement method is formulate by relating nodal displacements $\{q\}$ to the internal strain energy $U$. For a truss element with nodes $1$ and $2$, the change in length $\Delta L$ is the projection of the relative displacement onto the element's unit vector $\hat{n}$.

### 1. Kinematics
Let $\mathbf{u}_1$ and $\mathbf{u}_2$ ( both $\in \mathbb{R}^{n \times 1}$ )   be the displacement vectors of the two nodes. The axial deformation $\Delta L$ is:
$$\Delta L = (\mathbf{u}_2 - \mathbf{u}_1) \cdot \hat{n} = \begin{bmatrix} -\hat{n}^T & \hat{n}^T \end{bmatrix} \begin{Bmatrix} \mathbf{u}_1 \\ \mathbf{u}_2 \end{Bmatrix}$$

### 2. Strain Energy
The elastic strain energy for an axial member is $U = \frac{EA}{2L}(\Delta L)^2$. Expand $(\Delta L)^2$:
$$(\Delta L)^2 = \left( \begin{bmatrix} -\hat{n}^T & \hat{n}^T \end{bmatrix} \begin{Bmatrix} \mathbf{u}_1 \\ \mathbf{u}_2 \end{Bmatrix} \right)^T \left( \begin{bmatrix} -\hat{n}^T & \hat{n}^T \end{bmatrix} \begin{Bmatrix} \mathbf{u}_1 \\ \mathbf{u}_2 \end{Bmatrix} \right)$$

Expanding this results in the quadratic form:
$$(\Delta L)^2 = \begin{bmatrix} \mathbf{u}_1^T & \mathbf{u}_2^T \end{bmatrix} \underbrace{\begin{bmatrix} \hat{n}\hat{n}^T & -\hat{n}\hat{n}^T \\ -\hat{n}\hat{n}^T & \hat{n}\hat{n}^T \end{bmatrix}}_{\text{2*dim X 2*dim Matrix}} \begin{Bmatrix} \mathbf{u}_1 \\ \mathbf{u}_2 \end{Bmatrix}$$



Substituting this back into the energy equation $U = \frac{1}{2} \{q\}^T [k_e] \{q\}$:
$$U = \frac{1}{2} \{q\}^T \left( \frac{EA}{L} \begin{bmatrix} \hat{n}\hat{n}^T & -\hat{n}\hat{n}^T \\ -\hat{n}\hat{n}^T & \hat{n}\hat{n}^T \end{bmatrix} \right) \{q\}$$

Element stiffness matrix $[k_e]$ is defined :
$$[k_e] = \frac{EA}{L} \begin{bmatrix} \hat{n}\hat{n}^T & -\hat{n}\hat{n}^T \\ -\hat{n}\hat{n}^T & \hat{n}\hat{n}^T \end{bmatrix}$$

---

## Structure

The solver is designed with separate modules for linear algebra solvers and physics of the system.

### 1. Data Structure and Maths (`/components`)
-  **`element.h`**: Defines the structures (`Node`, `Element`, `Support`, `Load`). 
- **`solver.h`** : General linear algebra library. Still trying to optimise and add more features.
### 2. FEM (`/displacement`)
- **`assembly.h`**: Global assembly with `assembleGlobalStiffness()` and methods for finding individual element $[k_e]$ via ` elementStiffnessMatrix()` and local DOF mapping `elementDofMap()`.
- **`boundary.h`**: Implements boundary conditions using a row-column stripping method to maintain matrix symmetry.
- **`postprocess.h`**: Calculates Strain, Stress, Internal Force from the displacment obtained
    - `computeElementKinematics()` - Element's displacement $\delta_e$ along its axial unit vector $\hat{n}_e$.
    - `computeElementForce()` - Elements internal force $F_e = \frac{EA}{L}\Delta L$.
    - `computeAllElementResponses()`: The master post-processing function. It iterates through the entire element list and creates a vector of `ElementResponse` structs. For each element it calculates:
        1. **Strain** ($\epsilon$): $\frac{\Delta L}{L}$
        2. **Stress** ($\sigma$): $E \cdot \epsilon$
        3. **Force** ($F$): $A \cdot \sigma$

- **`displacement.h`** : Master FEM solver pipeline and solves the $Ku=F$ linear system with `conjugate_gradient()` method.
   

### 3. Verification (`/verification`)
Model integrity and safeguards against physically impossible or mathematically singular systems.

- **Indexing**:
    - `checkElementNodeValidity()`: Ensures every element uses a valid node index.
    - `checkDuplicateElements()`: Uses bit-shifting encoding on node ID's to detect and flag repeated element definitions.
- **Connectivity**:
  - `connectivityCheck()`: Stack DFS on a Static Adjacency List starting from the support nodes. If any node is unreachable, the system is flagged as disconnected, preventing a singular matrix error during the solve.
- **Determinacy & Stability**: 
    - `dof_constraint()`: It flags systems as under/over constrained/ Statically Determinate ($m + r < dim \cdot j$).
  
- `verifyModel()`: Performs all the above checks alove with the material properties ($E, A > 0$), loads dimensions and returns a `ValidationResults` struct with:
    - Overall result
    - Statically Determinate
    - Log message 
- Automatic detection of zero-length elements and singular matrix conditions inside the solver pipeline itself to avoid recomputations.

### 4. Static Equilibrium (`/static/statically_determined.h`)
Solves the truss system (if statically determinate, obtained via verification module) with Method of Joints to obtain internal forces.

- **`solve_truss()`**: Implements Method of Joints. Instead of a stiffness matrix, it constructs a Global Equilibrium Matrix $[A]$ based on the equlibiurm equations at nodes along each dimension.
- Reactions at supports are treated as additional unknowns in the system $[A]\{x\} = \{b\}$.
- This module is strictly for statically determinate systems. By comparing the forces from this module against the FEM results in `testing.cpp`, we can ensure the numerical integrity of the stiffness assembly.

---
### Verification (minimal)
`testing.cpp` includes tests for:
 - Mechanism detection
 - Disconnected Node 
 - Analytical check : single bar axial test
 - Cross method verification with static and displacmenet methods for a Warren truss.
---
**User-Defined ID's are decoupled from Internal memory offsets.
- No out of bounds errors occurs whe User IDs excees total node count (give `checkElementNodeValidity()` was also run).
- Frontend can use any naming convention
- In `postprocessing.h`, `node_i` and `node_j` use to access physical coordinates and `index_i`, `index_j` for specific row access in the solved displacement vector.
### Future plans
- Add input/output modules, and a frontend with WASM
- Transition to utilize sparse formats to handle larger systems more efficiently

