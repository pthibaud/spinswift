/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation

/// A class for managing a laser heating
/// 
/// Laser heating interactions and atomic effects. 
/// - Author: Pascal Thibaudeau
/// - Date: 30/07/2025
/// - Copyright: [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)
/// - Version: 0.1

class LaserExcitation : Codable {

    struct Temperatures : Codable {
        /// the electron temperature
        var electron : Double 
        /// the phonon or lattice temperature
        var phonon : Double 
        /// the spin temperature
        var spin : Double 

        init(electron : Double? = Double(), phonon: Double? = Double(), spin : Double? = Double()) {
            self.electron = electron!
            self.phonon = phonon!
            self.spin = spin!
        }

        /// addition of two temperatures
        static func + (a: Temperatures, b: Temperatures) -> Temperatures {
            return Temperatures(electron: (a.electron+b.electron), phonon: (a.phonon+b.phonon), spin: (a.spin+b.spin))
        }
        /// subtraction of two temperatures
        static func - (a: Temperatures, b: Temperatures) -> Temperatures {
            return Temperatures(electron: (a.electron-b.electron), phonon: (a.phonon-b.phonon), spin: (a.spin-b.spin))
        }
        /// addition and affectation of two temperatures
        static func += (a: inout Temperatures, b: Temperatures) {
            var c: Temperatures = Temperatures()
            c = a + b
            a = c
        }
        /// multiplication by a single scalar each of the temperatures
        static func * (a: Double, b: Temperatures) -> Temperatures {
            return Temperatures(electron: a*(b.electron), phonon: a*(b.phonon), spin: a*(b.spin))
        }
    }

    struct Pulse : Codable {
        /// nature of the laser pulse, either Gaussian or Triangle
        var form : String
        /// Laser fluence in J/cm²
        var fluence : Double
        var duration : Double
        var delay : Double

        init (form : String? = "Gaussian", fluence: Double? = Double(), duration: Double? = Double(), delay: Double? = Double()){
            self.form = form!
            self.fluence = fluence!
            self.duration = duration!
            self.delay = delay!
        }
    }

    struct TTM : Codable {

        typealias HeatCapacity = Temperatures

        struct Coupling : Codable {
            var electronPhonon : Double
            var electronSpin : Double
            var phononSpin : Double

            init(electronPhonon: Double? = Double(), electronSpin: Double? = Double(), phononSpin: Double? = Double()){
                self.electronPhonon = electronPhonon!
                self.electronSpin = electronSpin!
                self.phononSpin = phononSpin!
            }
        }
        var effectiveThickness : Double
        var initialTemperature : Double
        var damping : Double
        var heatCapacity : HeatCapacity?
        var coupling : Coupling?

        init (effectiveThickness: Double? = Double(), initialTemperature: Double? = Double(), damping: Double? = Double(), heatCapacity: HeatCapacity? = HeatCapacity(), coupling: Coupling? = Coupling()) {
            self.effectiveThickness = effectiveThickness!
            self.initialTemperature = initialTemperature!
            self.damping = damping!
            self.heatCapacity = heatCapacity!
            self.coupling = coupling!
        }
    }
    /// Output temperatures for the 3 or 2 temperatures model
    var temperatures : Temperatures
    /// Laser pulse parameters
    var pulse : Pulse
    /// Physical content of the 2-3 temperatures model
    var ttm : TTM
    /// Current time when the temperatures are evaluated
    var time : Double

    /**
        Initializes a new laser excitation with the provided parts and specifications.

        - Parameters:
            - time: Current time when the temperatures are evaluated
            - temperatures : structure of electron, phonon and spin temperatures
            - pulse : parameters structure of the laser pulse
            - ttm : structure of the physical content of the 2-3 temperatures model

        - Returns: A new laser excitation
    */
    
    init(time: Double? = Double(), temperatures: Temperatures? = Temperatures(), pulse: Pulse? = Pulse(), ttm : TTM? = TTM()) {
        self.time = time!
        self.temperatures = temperatures!
        self.pulse = pulse!
        self.ttm = ttm!
    }

    /// Compute the absorbed laser power in W/cm³
    func power(time: Double) -> Double {
        var power: Double = Double()

        switch self.pulse.form.lowercased() {
            case "gaussian":
                let Φ: Double = self.pulse.fluence
                let σ: Double = self.pulse.duration
                let δ: Double = self.pulse.delay
                let ζ: Double = self.ttm.effectiveThickness
                power = (Φ/(σ*ζ))*exp(-((time-δ)*(time-δ))/(0.36*σ*σ))
            default: break
        }
        return power
    }

    /// Compute the rate of variation of the temperatures in K/s
    private func rhs(time: Double, temperatures: Temperatures) -> Temperatures {
        let Cep: Double = self.ttm.coupling!.electronPhonon
        let Ces: Double = self.ttm.coupling!.electronSpin
        let Cps: Double = self.ttm.coupling!.phononSpin
        let γ: Double = self.ttm.heatCapacity!.electron // Ce=gamma*Te
        let Cp: Double = self.ttm.heatCapacity!.phonon
        let Cs: Double = self.ttm.heatCapacity!.spin
        let τ_ls: Double = self.ttm.damping
        let T_ref: Double = self.ttm.initialTemperature

        let Te : Double = temperatures.electron
        let Tp : Double = temperatures.phonon
        let Ts : Double = temperatures.spin

        var rate: Temperatures = LaserExcitation.Temperatures()

        var f0 : Double = power(time:time)
        f0 = f0/(γ*Te) // Laser power
        f0 -= (Cep/γ)*(1.0-(Tp/Te))
        f0 -= (Ces/γ)*(1.0-(Ts/Te))// f[0]=dTe/dt
        rate.electron = f0 - (1.0/(τ_ls))*(Te-T_ref) // Newton cooling
        rate.phonon = (Cep/Cp)*(Te-Tp)+(Cps/Cp)*(Ts-Tp) // f[1]=dTp/dt
        rate.spin = (Ces/Cs)*(Te-Ts)+(Cps/Cs)*(Tp-Ts) // f[2]=dTs/dt 
        return rate
    }

    /// Integration method of the 3 temperatures model for a Gaussian pulse
    func advanceTemperaturesGaussian(method:String, Δt : Double)  {
          switch method.lowercased() {
            /// Euler or Runge-Kutta of order 1
            case "euler", "rk1":
                self.temperatures += Δt*rhs(time:self.time,temperatures:self.temperatures)
            /// Runge-Kutta of order 2
            case "rk2":
                let k1: Temperatures = self.temperatures + 0.5*Δt*rhs(time:self.time,temperatures:self.temperatures)
                self.temperatures += Δt*rhs(time:self.time+0.5*Δt,temperatures:k1)
            /// Runge-Kutta of order 4
            case "rk4":
                let k1: Temperatures = rhs(time:self.time,temperatures:self.temperatures)
                let k2: Temperatures = rhs(time:self.time+0.5*Δt,temperatures:self.temperatures+0.5*Δt*k1)
                let k3: Temperatures = rhs(time:self.time+0.5*Δt,temperatures:self.temperatures+0.5*Δt*k2)
                let k4: Temperatures = rhs(time:self.time+Δt,temperatures:self.temperatures+Δt*k3)
                self.temperatures += (Δt/6)*(k1+2*k2+2*k3+k4)
            default: break
            }   
    }

    /// Estimation of the integration timestep needed for a given constant quality factor 0<Q<1
    func estimateTimestep(factor: Double) -> Double {
        let temp: LaserExcitation.Temperatures = rhs(time: self.time, temperatures: self.temperatures)
        let array: [Double] = [abs(temp.electron), abs(temp.phonon), abs(temp.spin)]
        return factor/(array.max()!)
    }

    func jsonify() throws -> String {
        let data: Data = try JSONEncoder().encode(self)
        let jsonString: String? = String(data:data, encoding:.utf8) 
        return jsonString!
    } 
}
