/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation

/// A class for managing atomic properties
///
/// An atomic system is a collection of atoms and relationships with them.
/// - Author: Pascal Thibaudeau
/// - Date: 14/04/2023
/// - Update author: Mouad Fattouhi
/// - Updated: 11/09/2024
/// - Version: 0.1
class Atom: Codable {
    /// Name of the atomic species
    var name: String
    /// Identifier of the atomic type
    var type: Int
    /// Atomic Landé factor in Bohr's magneton unit
    var g: Double
    /// Magnon energy
    var ℇ: Double
    /// Atomic cell volume
    var Vat: Double
    /// Exchange stifness at the critical temperature
    var Dref: Double
    /// van Hove singularity of magnon DOS
    var vanHove: Double
    /// Atomic Cartesian position
    var position: Vector3
    /// Atomic pulsation vector around its spin spins
    var ω: Vector3
    /// The first and second moments of dLLB
    var moments: Moments
    /** The first moment of the dLLB is the avergage of the spin <S> as a vector3
        The second moment of the dLLB is the statistical cumulant of the spin <S*S> which is a Matrix3
    **/

    /**
        Initializes a new atom with the provided parts and specifications.
    
        - Parameters:
            - name: The name of the atom
            - type: The type of the atom
            - position: The Cartesian position of the atom
            - spin: The spin direction of the atom as a unit vector
            - ω: The pulsation vector located on the atom
            - g: The Landé factor of the atom in Bohr's magneton unit
            - ℇ: The magnon energy in "mev"
            - moments: The first and second moments of dLLB <S_i> qnd <S_i*S_j> i,j=x,y,z
    
        - Returns: A new atom
    */
    init(
        name: String? = String(), type: Int? = Int(), position: Vector3? = Vector3(),
        ω: Vector3? = Vector3(), moments: Moments? = Moments(), g: Double? = Double(),
        ℇ: Double? = Double(), Vat: Double? = Double(), Dref: Double? = Double(),
        vanHove: Double? = Double()
    ) {
        //    if name.isEmpty {fatalError("A name has to be provided")}
        self.name = name!
        self.type = type!
        self.position = position!
        self.ω = ω!
        guard g! >= 0 else { fatalError("g factor must be positive!") }
        self.g = g!
        self.ℇ = ℇ!
        self.Vat = Vat!
        self.Dref = Dref!
        self.vanHove = vanHove!
        self.moments = moments!
    }

    struct Moments: Codable {
        /// This is the magnetic moment as a unit-vector where the components are $\langle spin_i(t) \rangle$
        var spin: Vector3
        /// Sigma is the second moment where the components are $\langle spin_i(t)\otimes spin_j(t)\rangle$
        var sigma: Matrix3

        init(spin: Vector3? = Vector3(), sigma: Matrix3? = Matrix3()) {
            self.spin = spin!
            self.sigma = sigma!
        }
        /// Addition of two moments
        static func + (a: Moments, b: Moments) -> Moments {
            return Moments(spin: (a.spin + b.spin), sigma: (a.sigma + b.sigma))
        }
        /// Accumumation of moments
        static func += (a: inout Moments, b: Moments) {
            var c: Moments = Moments()
            c = a + b
            a = c
        }
        /// multiplication of a moment by a scalar
        static func * (a: Double, b: Moments) -> Moments {
            return Moments(spin: a * (b.spin), sigma: a * (b.sigma))
        }
    }

    /*
        Compute the dynamics of the atomic spin first and second order cumulants using a set of 12 coupled equations.
    
        - Parameters:
          - method: The type of method used either [euler](euler), [RangeKutta](rk4 and rk8) or [Exponential integration](expI)
          - Δt: The timestep used (internal unit in ps).
    
        - Returns: The first and second order cumulants
    */

    private func computeThermalCoef(Temperature: Double? = Double(), thermostat: String) -> Double {
        let β: Double = (k_B.value * Temperature!)
        var coef: Double = 0.0
        switch thermostat.lowercased() {
        case "classical":
            coef = β
        case "quantum":
            coef = computeQFDRq(Temperature: Temperature!)
        default: coef = 0.0
        }

        return coef
    }

    private func rhs(
        moments: Moments, Temperature: Double? = Double(), alpha: Double? = Double(),
        thermostat: String
    ) -> Moments {
        let α: Double = alpha!
        let c: Double = 1 / (1 + (α * α))
        let n: Double = g * μ_B.value
        let D: Double =
            γ.value * (α / n)
            * computeThermalCoef(Temperature: Temperature!, thermostat: thermostat)
        let spin: Vector3 = moments.spin
        let Σ: Matrix3 = moments.sigma
        let One: Matrix3 = Matrix3(fill: "identity")
        var rhsdLLB: Moments = Atom.Moments()

        var a1: Matrix3 = (Σ.trace() * (ω ⊗ spin)) - (Σ.transpose() * (ω ⊗ spin))
        a1 += ((ω ⊗ spin) * Σ.transpose()) - ((ω ⊗ spin).trace() * Σ.transpose())
        a1 += ((ω ⊗ spin) - (ω ⊗ spin).transpose()) * Σ.transpose()
        a1 -= 2 * (((spin ⊗ spin).trace() * (ω ⊗ spin)) - ((spin ⊗ spin).transpose() * (ω ⊗ spin)))
        let m1: Matrix3 = (ω × Σ.transpose()) + (α * a1)
        let m2: Matrix3 = (2 * Σ.trace() * One) - (3 * (Σ + Σ.transpose()))

        rhsdLLB.spin =
            c * ((ω × spin) + ((α * Σ.trace() * ω) - (α * Σ.transpose() * ω)) - (2 * D * c * spin))
        rhsdLLB.sigma = c * (m1 + m1.transpose() + (D * c * m2))

        return rhsdLLB
    }

    /// advanceMoments integrates the new values of Moments in time
    /// - Parameters:
    ///   - method: Method of chose
    ///   - Δt: timestep
    ///   - T: temperature
    ///   - α: Gilbert damping value
    ///   - thermostat: nature of the thermostat
    ///
    func advanceMoments(
        method: String, Δt: Double, T: Double? = Double(), α: Double? = Double(),
        thermostat: String? = String()
    ) {
        switch method.lowercased() {
        case "euler":
            self.moments +=
                Δt * rhs(moments: self.moments, Temperature: T!, alpha: α!, thermostat: thermostat!)
        case "rk4":
            let k1: Moments = rhs(
                moments: self.moments, Temperature: T!, alpha: α!, thermostat: thermostat!)
            let k2: Moments = rhs(
                moments: self.moments + 0.5 * Δt * k1, Temperature: T!, alpha: α!,
                thermostat: thermostat!)
            let k3: Moments = rhs(
                moments: self.moments + 0.5 * Δt * k2, Temperature: T!, alpha: α!,
                thermostat: thermostat!)
            let k4: Moments = rhs(
                moments: self.moments + Δt * k3, Temperature: T!, alpha: α!, thermostat: thermostat!
            )
            self.moments += (Δt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        case "symplectic":
            print("Symplectic method is not implemented.")
        default: break
        }
    }

    private func computeQFDR(Temperature: Double) -> Double {
        let Tc: Double = 635
        let β: Double = (k_B.value * Temperature)
        var coef: Double = 0.0
        let D0: Double = 1 - (Temperature / Tc)
        let Ed: Double = ℇ * pow(D0, 1 / 3)
        let u: Double = (Ed / β)

        let x: [Double] = [
            -0.99726386, -0.98561151, -0.96476226, -0.93490608, -0.89632116, -0.84936761,
            -0.7944838, -0.73218212, -0.66304427, -0.58771576, -0.50689991, -0.42135128, -0.3318686,
            -0.23928736, -0.14447196, -0.04830767, 0.04830767, 0.14447196, 0.23928736, 0.3318686,
            0.42135128, 0.50689991, 0.58771576, 0.66304427, 0.73218212, 0.7944838, 0.84936761,
            0.89632116, 0.93490608, 0.96476226, 0.98561151, 0.99726386,
        ]

        let w: [Double] = [
            0.00701861, 0.01627439, 0.02539207, 0.03427386, 0.0428359, 0.05099806, 0.05868409,
            0.06582222, 0.07234579, 0.0781939, 0.08331192, 0.08765209, 0.09117388, 0.0938444,
            0.09563872, 0.09654009, 0.09654009, 0.09563872, 0.0938444, 0.09117388, 0.08765209,
            0.08331192, 0.0781939, 0.07234579, 0.06582222, 0.05868409, 0.05099806, 0.0428359,
            0.03427386, 0.02539207, 0.01627439, 0.00701861,
        ]

        if Temperature < Tc {
            for i: Int in 0..<32 {
                let xnew = 0.5 * u * (x[i] + 1)
                coef += w[i] * (pow(xnew, 1.5) / (exp(xnew) - 1))
            }
            coef *= 0.5 * u
            coef *= 1.5 * Ed * pow(u, -2.5)
        } else {
            coef = β
        }
        return coef
    }

    // m-DOS integrand
    private func mDOS(_ w: Double, T: Double, b: Double, D_T: Double) -> Double {

        let numerator = ℏ.value * w
        let denominator = exp(ℏ.value * w / (k_B.value * T)) - 1.0
        let sqrtInner = (1.0 - ((4.0 * ℏ.value * w * b) / D_T)).squareRoot()
        let term = ((1.0 / (2.0 * b)) * (1.0 - sqrtInner)).squareRoot() / sqrtInner

        return (numerator / denominator) * term
    }

    private func computeQFDRq(Temperature: Double) -> Double {
        let V0: Double = Vat
        var coef: Double = 0.0
        let D0: Double = ℇ
        let Dc: Double = Dref
        let b: Double = vanHove

        let D_T: Double = D0 * pow(moments.spin.norm(), 2) + Dc

        let pref = 1 * (ℏ.value * V0) / (4 * π * π * D_T)
        let w_c: Double = D_T / (4 * ℏ.value * b)  // cut-off frequency

        let x: [Double] = [
            -0.99726386, -0.98561151, -0.96476226, -0.93490608, -0.89632116, -0.84936761,
            -0.7944838, -0.73218212, -0.66304427, -0.58771576, -0.50689991, -0.42135128, -0.3318686,
            -0.23928736, -0.14447196, -0.04830767, 0.04830767, 0.14447196, 0.23928736, 0.3318686,
            0.42135128, 0.50689991, 0.58771576, 0.66304427, 0.73218212, 0.7944838, 0.84936761,
            0.89632116, 0.93490608, 0.96476226, 0.98561151, 0.99726386,
        ]

        let w: [Double] = [
            0.00701861, 0.01627439, 0.02539207, 0.03427386, 0.0428359, 0.05099806, 0.05868409,
            0.06582222, 0.07234579, 0.0781939, 0.08331192, 0.08765209, 0.09117388, 0.0938444,
            0.09563872, 0.09654009, 0.09654009, 0.09563872, 0.0938444, 0.09117388, 0.08765209,
            0.08331192, 0.0781939, 0.07234579, 0.06582222, 0.05868409, 0.05099806, 0.0428359,
            0.03427386, 0.02539207, 0.01627439, 0.00701861,
        ]
        for i: Int in 0..<w.count {
            let xnew = 0.5 * w_c * (x[i] + 1)
            coef += w[i] * mDOS(xnew, T: Temperature, b: b, D_T: D_T)
        }
        coef *= 0.5 * w_c * pref

        return coef
    }

    func jsonify() throws -> String {
        let data: Data = try JSONEncoder().encode(self)
        if let jsonString = String(data: data, encoding: .utf8) {
            return jsonString
        } else {
            throw NSError(
                domain: "EncodingError", code: 1,
                userInfo: [NSLocalizedDescriptionKey: "Failed to convert data to JSON string"])
        }
    }
}
