/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation
import CGSL

// atomistic code syntax
/*
var a: Atom = Atom(name: PeriodicTable().Z_label("Iron"), type:1, position: Vector3(1,1,1), spin: Vector3(direction:"+y"), g:2)
var b: Atom = Atom(name: PeriodicTable().Z_label("Iron"), type:1, position: Vector3(0,1,2), spin: Vector3(direction:"+y"), g:2)

var atoms: [Atom] = [a,b]

atoms.forEach {
    print(try! $0.jsonify())
    }
*/

// device syntax

var MTJ: Stack = Stack()
MTJ.append(Magnetization(name: "Fe", type: 1, position: Vector3(0,0,0), spin: Vector3(direction:"+x"), g:2))
MTJ.append(Magnetization(name: "O", type: 2, position: Vector3(0,0,1), spin: Vector3(0,0,0), g:0))
MTJ.append(Magnetization(name: "Fe", type: 1, position: Vector3(0,0,2), spin: Vector3(direction:"+x"), g:2))

var h: Interaction = Interaction([MTJ[0],MTJ[2]])
.zeemanField(Vector3(direction:"+z"), value: 0.1)
.uniaxialField(Vector3(direction:"+x"), value: 0.0)
//.exchangeField(typeI:1,typeJ:1,value:1.0,Rcut:3)
.dampening(0.1)

//print(try! h.jsonify())
let s: Integrate = Integrate(h)
s.expLs(Δt:0.1, by:"euler")
//print(try! h.atoms[0].jsonify())
//print(Analysis(h.atoms).temperature())

let pulse: LaserExcitation.Pulse = LaserExcitation.Pulse(form: "Gaussian", fluence: 10.0, duration: 60E-15, delay: 0)
let Cp: LaserExcitation.TTM.HeatCapacity = LaserExcitation.TTM.HeatCapacity(electron:6E3, phonon:2.2E6, spin:2.2E6)
let G: LaserExcitation.TTM.Coupling = LaserExcitation.TTM.Coupling(electronPhonon: 2.5E17)
let ttm: LaserExcitation.TTM = LaserExcitation.TTM(effectiveThickness: 1e-9, initialTemperature: 300, damping: 1E-12, heatCapacity: Cp, coupling: G)
let laser: LaserExcitation = LaserExcitation(temperatures: .init(electron: ttm.initialTemperature, phonon: ttm.initialTemperature, spin: ttm.initialTemperature), pulse: pulse,ttm: ttm)
//print(laser.power(time:laser.time))

for _ in 0...30000 {
    let timestep: Double = laser.estimateTimestep(factor:0.8)
    laser.advanceTemperaturesGaussian(method:"rk4",Δt:timestep)
    laser.time += timestep
    print(String(format:"%e %f %f %f",laser.time,laser.temperatures.electron,laser.temperatures.phonon,laser.temperatures.spin))
}
