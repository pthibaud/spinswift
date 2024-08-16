/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation
/// A class for integrating the dynamics of spins
/// 
/// The purpose of this class is to integrate a collection of particles. 
/// - Author: Pascal Thibaudeau
/// - Date: 03/10/2023
/// - Copyright: [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)
/// - Version: 0.1

class Integrate: Codable {

    var h: Interaction

    init(_ h: Interaction = Interaction()) {
        self.h = h
    }

    func evolve(_ iterations: Int, Δt: Double, by: String? = nil, file: String) {
        switch by?.lowercased(){
        case "euler"? :
        self.evolveEuler(iterations: iterations, Δt: Δt, fileName: file)
        case "expls"? :
        self.evolveExpLs(iterations: iterations, Δt: Δt, fileName: file)
        case "symplectic"? :
        self.evolveSymplectic(iterations: iterations, Δt: Δt, fileName: file)
        default: break
        }
    }

    private func evolveEuler(iterations: Int, Δt: Double, fileName: String) {
        var time: Double = 0.0
        var content: String = String()

        for _ in 1...iterations {
            content += String(time)

            for a: Atom in h.atoms {
               content += " "+String(a.spin.x)+" "+String(a.spin.y)+" "+String(a.spin.z)
            }

            DispatchQueue.concurrentPerform(iterations: h.atoms.count) { (index: Int) in
                h.atoms[index].advanceSpin(Δt: Δt, by: "euler")
            }

            content += "\n"
            h.update()
            time += Δt
        }
        //let home = FileManager.default.homeDirectoryForCurrentUser
        saveOnFile(data:content, fileName: fileName)
    }

    private func evolveExpLs(iterations: Int, Δt: Double, fileName: String) {

        var time: Double = 0.0
        var content: String = String()

        for _ in 1...iterations {
            content += String(time)

            for a: Atom in h.atoms {
               content += " "+String(a.spin.x)+" "+String(a.spin.y)+" "+String(a.spin.z)
               content += "\n"
            }

            expLs(Δt: Δt, by: "euler")

            time += Δt
        }
        //let home = FileManager.default.homeDirectoryForCurrentUser
        saveOnFile(data:content, fileName: fileName)

    }

    private func evolveSymplectic(iterations: Int, Δt: Double, fileName: String) {

        var time: Double = 0.0
        var content: String = String()

        for _ in 1...iterations {
            content += String(time)

            for a: Atom in h.atoms {
               content += " "+String(a.spin.x)+" "+String(a.spin.y)+" "+String(a.spin.z)
               content += "\n"
            }

            expLs(Δt: Δt, by: "symplectic")

            time += Δt
        }
        //let home = FileManager.default.homeDirectoryForCurrentUser
        saveOnFile(data:content, fileName: fileName)

    }

    private func expLs(Δt: Double, by: String) {
        let numberOfAtoms: Int = h.atoms.count
        for i: Int in 0...(numberOfAtoms-2) {
            h.atoms[i].advanceSpin(Δt: Δt/2, by: by)
            h.update()
        }
        h.atoms[numberOfAtoms-1].advanceSpin(Δt: Δt, by: by)
        h.update()
        for i: Int in 0...(numberOfAtoms-2) {
            h.atoms[numberOfAtoms-i-2].advanceSpin(Δt: Δt/2, by: by)
            h.update()
        }
    }

    func jsonify() throws -> String {
        let data: Data = try JSONEncoder().encode(self)
        let jsonString: String? = String(data:data, encoding:.utf8) 
        return jsonString!
    } 

}