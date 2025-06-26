/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation

/****** Simulate Nickel bulk FCC *******/

let a: Double = 0.35
let Jexp: Double = 17.2
let J_ij: Double = 0.79 * Jexp
let D_0: Double = J_ij * a * a
let V_0: Double = a * a * a
let Em: Double = D_0 * pow(((6 * π * π) / V_0), 2 / 3)
print(String(D_0))

//Define parameters

//Initialze dLLB moments
let mm = Atom.Moments(spin: Vector3(direction: "+x"), sigma: Matrix3(1, 0, 0, 0, 0, 0, 0, 0, 0))

//Initialze Atom
let initials = InitialParam(name: "Ni", type: 1, moments: mm, g: 2.02, ℇ: D_0)

//Inititialze simulation program
let inisimulation = SimulationProgram.Inputs(
    T_initial: 0, T_step: 50, T_final: 1800, time_step: 1e-2, stop: 10, α: 0.1,
    thermostat: "classical")

//Initialize primary crystal cell
let unitCellAtoms: [Atom] = [
    Atom(position: Vector3(0.0, 0.0, 0.0)),
    Atom(position: Vector3(0.0, 0.5, 0.5)),
    Atom(position: Vector3(0.5, 0.0, 0.5)),
    Atom(position: Vector3(0.5, 0.5, 0.0)),
]

// Define the supercell dimensions
let supercellDimensions = Vector3(3, 3, 3)

// Generate the cubic crystal structure
let latticeConstants = Vector3(a, a, a)
let crystalStructure = generateCrystalStructure(
    unitCellAtoms: unitCellAtoms, supercell: supercellDimensions,
    latticeConstants: latticeConstants,
    initialParam: initials)

//Define boundary conditions
let Boundaries = BoundaryConditions(
    boxSize: Vector3(
        latticeConstants.x * supercellDimensions.x,
        latticeConstants.y * supercellDimensions.y,
        latticeConstants.z * supercellDimensions.z),
    periodicity: "on")

//Set the desired iteraction (Here only exchange is used)
var h: Interaction = Interaction(crystalStructure)
    .exchangeField(typeI: 1, typeJ: 1, value: J_ij / ℏ.value, cutoffRadius: 0.25, BCs: Boundaries)
//.ZeemanField(Vector3(direction:"+z"), value: 1.5)

//Integrate the dLLB equation for each atom in thhe crystal
let sol: Integrate = Integrate(h)

//Simulation programe for which the integration will be done
var p = SimulationProgram(sol)
p.simulate(program: "curie_temperature", inputs: inisimulation)
