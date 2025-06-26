/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/

typealias Magnetization = Atom

typealias Stack = [Atom]

//Structure to define initial parameters for an atom
struct InitialParam: Codable {
    var name: String
    var type: Int
    var spin: Vector3
    var moments: Atom.Moments
    var position: Vector3
    var g: Double
    var ℇ: Double
    var Vat: Double
    var Dref: Double
    var vanHove: Double

    init(
        name: String, type: Int, spin: Vector3? = Vector3(), moments: Atom.Moments?,
        position: Vector3? = Vector3(), g: Double? = Double(), ℇ: Double? = Double(),
        Vat: Double? = Double(), Dref: Double? = Double(), vanHove: Double? = Double()
    ) {
        self.name = name
        self.type = type
        self.spin = spin!
        self.moments = moments!
        self.position = position!
        self.g = g!
        self.ℇ = ℇ!
        self.Vat = Vat!
        self.Dref = Dref!
        self.vanHove = vanHove!
    }
}

// Structure to define boundary conditions
struct BoundaryConditions: Codable {
    var boxSize: Vector3
    var periodicity: String

    init(boxSize: Vector3? = Vector3(), periodicity: String? = String()) {
        self.boxSize = boxSize!
        self.periodicity = periodicity!
    }
}

// Generate the crystal structures (change the initialization step)
func generateCrystalStructure(
    unitCellAtoms: [Atom], supercell: Vector3, latticeConstants: Vector3,
    initialParam: InitialParam
) -> [Atom] {
    var crystalStructure: [Atom] = []

    for i in 0..<Int(supercell.x) {
        for j in 0..<Int(supercell.y) {
            for k in 0..<Int(supercell.z) {
                let translationVector = Vector3(Double(i), Double(j), Double(k))
                unitCellAtoms.forEach { atom in
                    let newPosition = Vector3(
                        latticeConstants.x * (atom.position.x + translationVector.x),
                        latticeConstants.y * (atom.position.y + translationVector.y),
                        latticeConstants.z * (atom.position.z + translationVector.z)
                    )
                    let newAtom = createAtom(position: newPosition, initialParam: initialParam)
                    crystalStructure.append(newAtom)
                }
            }
        }
    }

    return crystalStructure
}

private func createAtom(position: Vector3, initialParam: InitialParam) -> Atom {
    let newAtom = Atom(position: position)
    newAtom.name = initialParam.name
    newAtom.type = initialParam.type
    newAtom.moments = initialParam.moments
    newAtom.g = initialParam.g
    newAtom.ℇ = initialParam.ℇ
    newAtom.Vat = initialParam.Vat
    newAtom.Dref = initialParam.Dref
    newAtom.vanHove = initialParam.vanHove
    return newAtom
}

func substituteRandomAtoms(structure: [Atom], InitParam: InitialParam, Percentage: Double) -> [Atom]
{
    let Alloy: [Atom] = structure
    let N: Double = Percentage / 100  // round it is better
    let atomsSubstituted = N * Double(Alloy.count)
    let RandomAtoms = Array(0...Alloy.count - 1).shuffled()

    for i in 0..<Int(atomsSubstituted) {
        let atomicIndex = RandomAtoms[i]
        Alloy[atomicIndex].name = InitParam.name
        Alloy[atomicIndex].type = InitParam.type
        Alloy[atomicIndex].moments = InitParam.moments
        Alloy[atomicIndex].g = InitParam.g
        Alloy[atomicIndex].ℇ = InitParam.ℇ
        Alloy[atomicIndex].Vat = InitParam.Vat
        Alloy[atomicIndex].Dref = InitParam.Dref
        Alloy[atomicIndex].vanHove = InitParam.vanHove
    }
    return Alloy
}

func substituteSpecificAtoms(
    structure: [Atom], InitParam: InitialParam, unitCellAtoms: [Atom],
    supercellSize: Vector3
) -> [Atom] {
    // Create a deep copy of the original structure
    let Alloy = structure

    // Identify the reference atom in the unit cell (e.g., the one at (0,0,0))
    guard let referenceAtom = unitCellAtoms.first else {
        print("Error: Unit cell is empty.")
        return Alloy
    }

    // Function to check if an atom is a translation of the reference atom
    func isTranslation(of reference: Atom, atom: Atom, supercellSize: Vector3) -> Bool {

        return (atom.position.x - reference.position.x).truncatingRemainder(
            dividingBy: supercellSize.x)
            == 0
            && (atom.position.y - reference.position.y).truncatingRemainder(
                dividingBy: supercellSize.y)
                == 0
            && (atom.position.z - reference.position.z).truncatingRemainder(
                dividingBy: supercellSize.z)
                == 0
    }

    // Perform substitution only on atoms that match the reference pattern
    for atom in Alloy {
        if isTranslation(of: referenceAtom, atom: atom, supercellSize: supercellSize) {
            atom.name = InitParam.name
            atom.type = InitParam.type
            atom.moments = InitParam.moments
            atom.g = InitParam.g
            atom.ℇ = InitParam.ℇ
            atom.Vat = InitParam.Vat
            atom.Dref = InitParam.Dref
            atom.vanHove = InitParam.vanHove
        }
    }

    return Alloy
}

//Calculate the inter-atomic distances with or without periodic boundary conditions
func computeDistance(boundaryConditions: BoundaryConditions, atom1: Atom, atom2: Atom) -> Double {
    let boxSize = boundaryConditions.boxSize
    let periodicity = boundaryConditions.periodicity.lowercased()

    switch periodicity {
    case "off":
        return distance(atom1.position, atom2.position)
    case "on":
        var xij = atom1.position - atom2.position
        xij.x -= boxSize.x * (xij.x / boxSize.x).rounded(.toNearestOrAwayFromZero)
        xij.y -= boxSize.y * (xij.y / boxSize.y).rounded(.toNearestOrAwayFromZero)
        xij.z -= boxSize.z * (xij.z / boxSize.z).rounded(.toNearestOrAwayFromZero)
        return (xij ° xij).squareRoot()
    default:
        return 0.0
    }
}
