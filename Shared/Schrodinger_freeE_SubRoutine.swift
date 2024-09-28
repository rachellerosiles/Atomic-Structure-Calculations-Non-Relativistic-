//
//  Schrodinger_freeE_SubRoutine.swift
//  Atomic Structure Calculations (Non-Relativistic) (iOS)
//
//  Created by Rachelle Rosiles on 8/19/24.
//

import Foundation

class WaveFunctionCalculator {
    // Parameters
    var z_value: Double
    var trial_energy: Double
    var number_blocks: Double
    var l_number: Double
    var number_points: Double
    var initial_delta_x: Double
    var mesh_scalar: Double
    var principal_quant_number: Double
    var input_pot_list: [Double]
    var thresh_criterion: Double
    var max_thresh_iterations: Int
    var left_energy_scalar: Double
    var right_energy_scalar: Double
    
    init(z_value: Double, trial_energy: Double, number_blocks: Double, l_number: Double,
         number_points: Double, initial_delta_x: Double, mesh_scalar: Double,
         principal_quant_number: Double, input_pot_list: [Double], thresh_criterion: Double,
         max_thresh_iterations: Int, left_energy_scalar: Double, right_energy_scalar: Double) {
        self.z_value = z_value
        self.trial_energy = trial_energy
        self.number_blocks = number_blocks
        self.l_number = l_number
        self.number_points = number_points
        self.initial_delta_x = initial_delta_x
        self.mesh_scalar = mesh_scalar
        self.principal_quant_number = principal_quant_number
        self.input_pot_list = input_pot_list
        self.thresh_criterion = thresh_criterion
        self.max_thresh_iterations = max_thresh_iterations
        self.left_energy_scalar = left_energy_scalar
        self.right_energy_scalar = right_energy_scalar
    }
    
    //toy V example
    func V(r: Double, customPotential: ((Double) -> Double)? = nil) -> Double {
        if r == 0 {
            return Double.infinity // Avoid singularity
        }
        if let potential = customPotential {
            return potential(r)
        }
        return -z_value / r // Default Coulomb potential
    }
    
    func derivatives(r: Double, y: [Double]) -> [Double] {
        let y1 = y[0] // P(r)
        let y2 = y[1] // P'(r)
        
        let dY1 = y2
        let dY2 = -1 * (V(r: r) + trial_energy - (l_number * (l_number + 1)) / (r * r) * y1)
        
        return [dY1, dY2]
    }
    
    //Runge-Kutta 4th order method
    func rungeKutta(h: Double, r0: Double, y0: [Double], rEnd: Double) -> [(Double, [Double])] {
        var r = r0
        var y = y0
        var results: [(Double, [Double])] = []

        while r <= rEnd {
            results.append((r, y))

            let k1 = derivatives(r: r, y: y)

            let y2Temp = zip(y, k1).map { $0 + $1 * h / 2.0 }
            let k2 = derivatives(r: r + h / 2.0, y: y2Temp)

            let y3Temp = zip(y, k2).map { $0 + $1 * h / 2.0 }
            let k3 = derivatives(r: r + h / 2.0, y: y3Temp)

            let y4Temp = zip(y, k3).map { $0 + $1 * h }
            let k4 = derivatives(r: r + h, y: y4Temp)

            // Update y using the RK4
            let yNew = zip(y, zip(k1, zip(k2, zip(k3, k4))).map { (k1, k2, k3, k4) in
                return k1 + 2 * k2 + 2 * k3 + k4
            }).map { (yValue, kSum) in
                return yValue + (kSum * h / 6.0)
            }

            y = yNew
            r += h
        }

        return results
    }
}
