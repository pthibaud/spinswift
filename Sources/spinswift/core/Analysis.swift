/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation

/// A class for getting physical properties from atomic spin states
/// 
/// - Author: Pascal Thibaudeau
/// - Date: 14/04/2024
/// - Copyright: [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)
/// - Version: 0.1
class Analysis {

    var atoms : [Atom] 

    init(_ atoms: [Atom]? = [Atom]()) {
        self.atoms = atoms!
    }

    /// Return the total energy of a magnetic system
    func GetEnergy() -> Double {
        var e : Double = 0
        atoms.forEach {
            e += ($0.ω°$0.spin)
        }
        return e
    }

    /// Return the average magnetization vector of a magnetic system
    func GetMagnetization() -> Vector3 {
        var m : Vector3 = Vector3()
        var g : Double = 0
        atoms.forEach {
            m += (($0.g)*($0.spin))
            g += ($0.g)
        }

        guard g != 0 else {return Vector3(0,0,0)}
        return (1.0/g)*m
    }

    /// Return the total torque of a magnetic system
    func GetTorque() -> Vector3 {
        var t: Vector3 = Vector3()
        atoms.forEach {
            t += ($0.ω)×($0.spin)
        }
        return t
    }

    /// Return the spin microcanonic temperature 
    func GetTemperature(coefficient:Double? = 2.0) -> Double {
        let e: Double = self.GetEnergy()
        let t: Vector3 = self.GetTorque()
        let t2: Double = t°t
        let T: Double = (t2*(ℏ.value))/(e*coefficient!*(k_B.value))
        return T
    }
}
