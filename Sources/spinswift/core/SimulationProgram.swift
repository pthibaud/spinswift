/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation
/// This class for managing The simualtion programs of spinswift code
/// All simualtion programs are created and managed in this class, the list of available simulation programs is bellow
/// - Author: Mouad Fattouhi 
/// - Date: 02/09/2024
/// - Version: 0.1

/* List of simulation programs:

   1) Curie_Temperature: This program allows for simulating multispin dynamics with sweeped temprerature within a range provided by the user in order to determine the Curie curve. 

   2) Optical_Pulse: This Program simulates the effects of an applied laser pulse on magnetization by conecting a 2 or 3 temperatures model coupled to the
      magnetization dynamics equations.

   3) time_dynamics: This simualtion program allows for simualting a single spin dynamics in presence of user defined magnetic interaction.

   4) ParamgneticSpins: This program allows for simulating a the dynamics of a multispin system in abscence of exhcange interactions. 
*/

class SimulationProgram : Codable {

    var I:Integrate

    init(_ I: Integrate = Integrate()) {
        self.I = I
    }

    struct Inputs : Codable {
        var T_initial: Double
        var T_step: Double
        var T_final: Double
        var stop: Double
        var time_step: Double
        var α : Double
        var thermostat : String
        
        init(T_initial: Double? = Double(), T_step: Double? = Double(), T_final: Double? = Double(), time_step: Double? = Double(), stop: Double? = Double(), α: Double? = Double(), thermostat: String? = String()){
            self.T_initial = T_initial!
            self.T_step = T_step! 
            self.T_final = T_final!
            self.time_step = time_step!
            self.stop = stop!
            self.α = α!
            self.thermostat = thermostat!    
        }
    }

    func simulate(program: String? = nil, inputs: Inputs) {
        guard let program = program else {
            print("Program not found. Choose one of the available programs:\ncurie_temperature\noptical_pulse\nmacrospin\nparamagneticspins")
            return
        }

        switch program.lowercased() {
        case "curie_temperature":
            curieTemperature(initialize: inputs)
        case "optical_pulse":
            opticalPulse(inputs: inputs)
        case "time_dynamics":
            timeDynamics(initialize: inputs)
        case "paramagneticspins":
            break
        }
    }

    //Curie temperature program

    private func curieTemperature(Initialize: Inputs) {
        var T: Double = Initialize.T_initial
        var content1: String = String()
        var content2: String = String()
        let stop: Double = Initialize.stop; let Δt : Double = Initialize.time_step; let T_final: Double = Initialize.T_final
        let dT: Double = Initialize.T_step; let α: Double = Initialize.α; let thermostat: String = Initialize.thermostat 
        var sl1: [Atom] = []
        var sl2: [Atom] = []

       for i in I.h.atoms {
                if (i.type == 1) {sl1.append(i)}
                else {sl2.append(i)}
            }

        while (T < T_final) {
            let Fn : String = "CT_T_" + String(format: "%.0f", T)
            I.evolveRK45(stop: stop, Δt: Δt, fileName: Fn, Temp: T, alpha: α, thermostat: thermostat)
            let m : Vector3 = Analysis(I.h.atoms).getMagnetization()
            let m1 : Vector3 = Analysis(sl1).getMagnetization()
            let m2 : Vector3 = Analysis(sl2).getMagnetization()
            let mnorm : Double = Analysis(I.h.atoms).getMagnetizationLength()
            let mnorm1 : Double = Analysis(sl1).getMagnetizationLength()
            let mnorm2 : Double = Analysis(sl2).getMagnetizationLength()
            let χ : Matrix3 = Analysis(I.h.atoms).getSusceptibility()
            content1 += String(T)+"\t"+String(m.x)+"\t"+String(m.y)+"\t"+String(m.z)+"\t"+String(mnorm)+"\t"
            content1 += String(χ.xx)+"\t"+String(χ.yy)+"\t"+String(χ.zz)+"\t"
            content1 += String(χ.xy)+"\t"+String(χ.yz)+"\t"+String(χ.zx)+"\n"

            content2 += String(T)+"\t"+String(sl1[1].g*mnorm1)+"\t"+String(sl2[1].g*mnorm2)+"\t"+String(mnorm)+"\t"
            content2 += String(m1.z)+"\t"+String(m2.z)+"\t"+String(m.z)+"\n"
            T+=dT
        } 
        //k_B.value*
        SaveOnFile(data:content1, fileName: "Output_CurieTemp_MvsTfile_Gd_tst")
        SaveOnFile(data:content2, fileName: "Output_CurieTemp_MvsTfile_FeGd333")
    }

    //Laser Pulse programe

   private func opticalPulse(inputs: Inputs) {
        let pulse = LaserExcitation.Pulse(Form: "Gaussian", Fluence: 32.5, Duration: 60E-15, Delay: 5e-12)
        let Cp = LaserExcitation.TTM.HeatCapacity(Electron:7E3, Phonon:3e6)
        let G = LaserExcitation.TTM.Coupling(ElectronPhonon: 60e17)
        let ttm = LaserExcitation.TTM(EffectiveThickness:15E-9, InitialTemperature: 82, Damping: 5E-12, HeatCapacity: Cp, Coupling: G)
        let laser = LaserExcitation(temperatures: .init(Electron:ttm.InitialTemperature, Phonon:ttm.InitialTemperature, Spin:ttm.InitialTemperature), pulse:pulse,ttm:ttm)
        var content: String = String()
        let Δt : Double = inputs.time_step
        var sl1: [Atom] = []
        var sl2: [Atom] = []

        for i in I.h.atoms {
                if (i.type == 1) {sl1.append(i)}
                else {sl2.append(i)}
            }

        while (laser.CurrentTime < inputs.stop*1E-12) {
            laser.AdvanceTemperaturesGaussian(method:"euler",Δt: Δt*1E-12)
            laser.CurrentTime += Δt*1E-12
            for a in I.h.atoms {
                a.advanceMoments(method: "rk4", Δt: Δt, T: laser.temperatures.Electron, α: inputs.α, thermostat: inputs.thermostat)
            } 

            //laser.temperatures.Electron
            let m : Vector3 = Analysis(I.h.atoms).getMagnetization()
            let m1 : Vector3 = Analysis(sl1).getMagnetization()
            let m2 : Vector3 = Analysis(sl2).getMagnetization()
            let mnorm : Double = Analysis(I.h.atoms).getMagnetizationLength()
            let mnorm1 : Double = Analysis(sl1).getMagnetizationLength()
            let mnorm2 : Double = Analysis(sl2).getMagnetizationLength()
            content += String(laser.CurrentTime)+"\t"+String(m1.z)+"\t"+String(m2.z)+"\t"+String(m.z)+"\t"
            content += String(laser.ComputeInstantPower(time:laser.CurrentTime))+"\t"+String(laser.temperatures.Electron)+"\t"+String(laser.temperatures.Phonon)+"\t"+String(laser.temperatures.Spin)+"\n"

            self.I.h.Update() 
        }
        SaveOnFile(data:content, fileName: "AOS_FeGd_25pct")
    }

    private func timeDynamics(Initialize: Inputs) { 
        let T: Double = Initialize.T_initial
        let Fn: String = "Dy_T_" + String(format: "%.0f", T) + "a_1e-2"
        let stop: Double = Initialize.stop; let Δt : Double = Initialize.time_step; let T_final: Double = Initialize.T_final
        let dT: Double = Initialize.T_step; let α: Double = Initialize.α; let thermostat: String = Initialize.thermostat

        I.evolveRK45(stop: stop, Δt: Δt, fileName: Fn, Temp: T, alpha: α, thermostat: thermostat)
    }
}

