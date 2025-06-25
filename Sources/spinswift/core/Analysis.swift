/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/

class Analysis {

    var atoms: [Atom]

    init(_ atoms: [Atom]? = nil) {
        self.atoms = atoms ?? []
    }

    func getEnergy() -> Double {
        return atoms.reduce(0) { $0 + ($1.ω ° $1.moments.spin) }
    }

    func getMagnetization() -> Vector3 {
        let (m, g) = atoms.reduce((Vector3.zero, 0)) { (acc, atom) in
            return (acc.0 + (atom.g * atom.moments.spin), acc.1 + atom.g)
        }
        guard g != 0 else { return .zero }
        return (1.0 / g) * m
    }

    func getMagnetizationLength() -> Double {
        let totalG = atoms.reduce(0) { $0 + $1.g }
        guard totalG != 0 else { return 0 }

        let mnorm = atoms.reduce(0) { $0 + $1.g * $1.moments.spin.norm() }
        return mnorm / totalG
    }

    func getTorque() -> Vector3 {
        return atoms.reduce(Vector3.zero) { $0 + ($1.ω × ($1.moments.spin)) }
    }

    func getTemperature(coefficient: Double? = nil) -> Double {
        let e = self.getEnergy()
        let t = self.getTorque()
        let t2 = t ° t
        let T = (t2 * ℏ.value) / (e * (coefficient ?? 2.0) * k_B.value)
        return T
    }

    func getInstantEnergy() -> Double {
        return atoms.reduce(0) { $0 + ($1.ω * $1.moments.spin) }
    }

    func getSusceptibility() -> Matrix3 {
        let N = Double(atoms.count)
        guard N != 0 else { return Matrix3() }

        let χ = atoms.reduce(Matrix3()) { result, atom in
            let A = atom.moments.spin ⊗ atom.moments.spin
            return result + atom.moments.sigma - A
        }

        return (1/N) * χ
    }

    func getCumulant() -> Matrix3 {
        let N = Double(atoms.count)
        guard N != 0 else { return Matrix3() }

        let Σ = atoms.reduce(Matrix3()) { result, atom in
            return result + atom.moments.sigma
        }

        return (1/N) * Σ
    }
}
