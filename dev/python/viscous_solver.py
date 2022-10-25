from dataclasses import dataclass
from math import fabs, sqrt, tanh
from typing import Iterable
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#----------------------------------------------#
# Constants
#----------------------------------------------#
layers = 300

eta = np.linspace(0, 1, num=layers)

f0 = 6 * np.power(eta, 2) - 8 * np.power(eta, 3) + 3 * np.power(eta, 4)
f1 = eta - 3 * np.power(eta, 2) + 3 * np.power(eta, 3) - np.power(eta, 4)
f2 = (eta - 4 * np.power(eta, 2) + 6 * np.power(eta, 3) - 4 * np.power(eta, 4) + np.power(eta, 5)) * np.power(1 - eta, 2)
f3 = (np.power(eta, 2) - 3 * np.power(eta, 3) + 3 * np.power(eta, 4) - np.power(eta, 5)) * np.power(1 - eta, 2)

df0_deta = 12 * eta - 24 * np.power(eta, 2) + 12 * np.power(eta, 3)
df1_deta = 1 - 6 * eta + 9 * np.power(eta, 2) - 4 * np.power(eta, 3)
df2_deta = (1 - 6 * eta + 9 * np.power(eta, 2) - 4 * np.power(eta, 3)) * np.power(1 - eta, 2) + 2 * (eta - 4 * np.power(eta, 2) + 6 * np.power(eta, 3) - 4 * np.power(eta, 4) + np.power(eta, 5)) * (1 - eta)
df3_deta = (2 * eta - 9 * np.power(eta, 2) + 12 * np.power(eta, 3) - 5 * np.power(eta, 4)) * np.power(1 - eta, 2) + 2 * (np.power(eta, 2) - 3 * np.power(eta, 3) + 3 * np.power(eta, 4) - np.power(eta, 5)) * (1 - eta)

#----------------------------------------------#
# Classes
#----------------------------------------------#
@dataclass
class Freestream:
    velocity: float = 10.0
    density: float = 1.225
    viscosity: float = 1e-5
    mach: float = 0.01

    def to_list(self) -> Iterable[float]:
        return (self.velocity, self.density, self.viscosity, self.mach)

@dataclass
class Divergentes:
    div_phi_11_12: float = 0.001
    div_phi_21_22: float = 0.001
    div_delta_1_2_ast: float = 0.001
    div_phi_1_2_ast: float = 0.001
    div_theta_1_2_o: float = 0.001
    div_Ktau_xx_xy: float = 0.001
    div_Ktau_yx_yy: float = 0.001

    def to_list(self) -> Iterable[float]:
        return (self.div_phi_11_12, self.div_phi_21_22, self.div_delta_1_2_ast, self.div_phi_1_2_ast, self.div_theta_1_2_o, self.div_Ktau_xx_xy, self.div_Ktau_yx_yy)

@dataclass
class Gradients:
    grad_u2_x: float = 0.001
    grad_u2_y: float = 0.001
    grad_phi_x: float = 0.001
    grad_phi_y: float = 0.001

    def to_list(self) -> Iterable[float]:
        return (self.grad_u2_x, self.grad_u2_y, self.grad_phi_x, self.grad_phi_y)

@dataclass
class Parameters:
    delta: float = 1e-5
    A: float = 1e-5
    B: float = 0.0
    Psi: float = 0.0
    Ctau1: float = 1e-15
    Ctau2: float = 1e-15
    eps: float = 1e-12

    delta_norm: float = 1.0
    A_norm: float = 1.0
    B_norm: float = 1.0
    Psi_norm: float = 1.0
    Ctau1_norm: float = 1.0
    Ctau2_norm: float = 1.0

    def add(self, values: Iterable[float]) -> None:
        self.delta = self.delta + values[0]
        self.A = self.A + values[1]
        self.B = self.B + values[2]
        self.Psi = self.Psi + values[3]
        self.Ctau1 = self.Ctau1 + values[4]
        self.Ctau2 = self.Ctau2 + values[5]
        return

    def to_list_norm(self) -> Iterable[float]:
        return (self.delta_norm, self.A_norm, self.B_norm, self.Psi_norm, self.Ctau1_norm, self.Ctau2_norm)
    
    def to_list(self) -> Iterable[float]:
        return (self.delta_norm * self.delta, self.A_norm * self.A, self.B_norm * self.B, self.Psi_norm * self.Psi, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * self.Ctau2)
    
    def to_list_eps(self, name: str) -> Iterable[float]:
        if name == 'delta':
            return (self.delta_norm * (self.delta) + self.eps, self.A_norm * self.A, self.B_norm * self.B, self.Psi_norm * self.Psi, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * self.Ctau2)
        elif name == 'A':
            return (self.delta_norm * self.delta, self.A_norm * (self.A) + self.eps, self.B_norm * self.B, self.Psi_norm * self.Psi, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * self.Ctau2)
        elif name == 'B':
            return (self.delta_norm * self.delta, self.A_norm * self.A, self.B_norm * (self.B) + self.eps, self.Psi_norm * self.Psi, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * self.Ctau2)
        elif name == 'Psi':
            return (self.delta_norm * self.delta, self.A_norm * self.A, self.B_norm * self.B, self.Psi_norm * (self.Psi) + self.eps, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * self.Ctau2)
        elif name == 'Ctau1':
            return (self.delta_norm * self.delta, self.A_norm * self.A, self.B_norm * self.B, self.Psi_norm * self.Psi, self.Ctau1_norm * (self.Ctau1) + self.eps, self.Ctau2_norm * self.Ctau2)
        elif name == 'Ctau2':
            return (self.delta_norm * self.delta, self.A_norm * self.A, self.B_norm * self.B, self.Psi_norm * self.Psi, self.Ctau1_norm * self.Ctau1, self.Ctau2_norm * (self.Ctau2) + self.eps)
        else:
            return ()
    
    def convert_to_real(self, x: Iterable[float]) -> Iterable[float]:
        return [x[0] / self.delta_norm, x[1] / self.A_norm, x[2] / self.B_norm, x[3] / self.Psi_norm, x[4] / self.Ctau1_norm, x[5] / self.Ctau2_norm]

#----------------------------------------------#
# Profiles
#----------------------------------------------#
def profiles_func(x: Iterable[float], air: Iterable[float]) -> Iterable[float]:

    delta, A, B, Psi, _, _ = x
    velocity, density, viscosity, mach = air

    # Velocities
    U = f0 + A * (1 - 0.6 * (A - 3) * eta * eta * eta) * f1
    W = B * f2 + Psi * f3

    # Velocity derivatives
    dU_deta = df0_deta + A * (1 - 0.6 * (A - 3) * eta * eta * eta) * df1_deta - A * 1.8 * (A - 3) * np.power(eta, 2) * f1
    dW_deta = B * df2_deta + Psi * df3_deta

    # Density
    R = 1 / (1 + 0.2 * mach * mach * (1 - U * U - W * W))

    # Viscosity
    MU = np.power(1 / R, 1.5) * (1 + 1 / R) / (1 / R + 1 / R)

    # Reynolds (delta)
    reynolds = velocity * density * delta / viscosity

    # Shear stress
    S = (1 / (reynolds + 1e-12)) * MU * dU_deta
    T = (1 / (reynolds + 1e-12)) * MU * dW_deta

    return [eta, U, W, dU_deta, dW_deta, R, S, T]

#----------------------------------------------#
# Integral thickness
#----------------------------------------------#
def integral_thickness_func(x: Iterable[float], profiles: Iterable[float], air: Iterable[float]) -> Iterable[float]:

    delta, A, B, Psi, Ctau1, Ctau2 = x
    eta, U, W, dU_deta, dW_deta, R, S, T = profiles
    velocity, density, viscosity, mach = air

    delta_1_ast = delta * simpson(1 - R * U, eta)
    delta_2_ast = delta * simpson(- R * W, eta)

    phi_11 = delta * simpson(1 - R * U * U, eta)
    phi_12 = delta * simpson(- R * U * W, eta)
    phi_21 = phi_12
    phi_22 = delta * simpson(- R * W * W, eta)

    phi_1_ast = delta * simpson(1 - R * U * (U * U + W * W), eta)
    phi_2_ast = delta * simpson(- R * U * (U * U + W * W), eta)

    delta_1_line = delta * simpson(1 - U, eta)
    delta_2_line = delta * simpson(- W, eta)

    delta_q = phi_11 + phi_22

    delta_q_o = delta * Psi * simpson(- R * (U * U + W * W), eta)

    theta_1_o = delta * Psi * simpson(- R * U * (U * U + W * W), eta)
    theta_2_o = delta * Psi * simpson(- R * W * (U * U + W * W), eta)

    delta_1_o = delta * Psi * simpson(- U, eta)
    delta_2_o = delta * Psi * simpson(- W, eta)

    S0 = S[0]
    T0 = T[0]

    C_D = simpson(S * dU_deta + T * dW_deta, eta)
    C_D_x = simpson(S * dW_deta - T * dU_deta, eta)
    C_D_o = Psi * simpson(S * dU_deta + T * dW_deta, eta)

    theta_11 = phi_11 - delta_1_ast

    H_k = (delta_1_ast / (phi_11 - delta_1_ast)  - 0.29 * mach * mach) / (1 + 0.113 * mach * mach)

    if (H_k < 1.8):

        Sx = 0.0
            
    else:

        Re_theta = velocity * density * abs(phi_11 - delta_1_ast) / viscosity
        func1 = 0.01 * sqrt(pow(2.4 * H_k - 3.7 + 2.5 * tanh(1.5 * H_k - 4.65), 2) + 0.25)
        func2 = pow(10, (1.415 / (H_k - 1) - 0.489) * tanh(20 / (H_k - 1) - 12.9) + 3.295 / (H_k - 1) + 0.44)
        fn = func1 * (Re_theta - func2)

        if (fn < 0):
            Sx = 0.0
        else:
            Sx = sqrt(Ctau1 * Ctau1 + Ctau2 * Ctau2) * fn / ((phi_11 - delta_1_ast) * density * velocity * velocity)
        
    H_k = (delta_2_ast / (phi_11 - delta_1_ast)  - 0.29 * mach * mach) / (1 + 0.113 * mach * mach)

    if (H_k < 1.8):
            
        Sy = 0.0
            
    else:

        Re_theta = velocity * density * abs(phi_22 - delta_2_ast) / viscosity
        func1 = 0.01 * sqrt(pow(2.4 * H_k - 3.7 + 2.5 * tanh(1.5 * H_k - 4.65), 2) + 0.25)
        func2 = pow(10, (1.415 / (H_k - 1) - 0.489) * tanh(20 / (H_k - 1) - 12.9) + 3.295 / (H_k - 1) + 0.44)
        fn = func1 * (Re_theta - func2)

        if (fn < 0):
            Sy = 0.0
        else:
            Sy = sqrt(Ctau1 * Ctau1 + Ctau2 * Ctau2) * fn / ((phi_11 - delta_1_ast) * density * velocity * velocity)
            

    return [delta_1_ast, delta_2_ast, phi_11, phi_12, phi_21, phi_22, phi_1_ast, phi_2_ast, delta_1_line, delta_2_line, delta_q, delta_q_o, theta_1_o, theta_2_o, delta_1_o, delta_2_o, S0, T0, C_D, C_D_x, C_D_o, theta_11, Sx, Sy]

#----------------------------------------------#
# Equations
#----------------------------------------------#
def equations_func(x: Iterable[float], divs: Iterable[float], grads: Iterable[float], air: Iterable[float]) -> Iterable[float]:

    # Profiles
    profiles = profiles_func(x, air)

    # Integral thickness
    params = integral_thickness_func(x, profiles, air)

    # Assign parameters
    _, _, _, _, _, _, phi_1_ast, phi_2_ast, delta_1_line, delta_2_line, _, _, _, _, delta_1_o, delta_2_o, S0, T0, C_D, C_D_x, C_D_o, _, Sx, Sy = params
    div_phi_11_12, div_phi_21_22, div_delta_1_2_ast, div_phi_1_2_ast, div_theta_1_2_o, div_Ktau_xx_xy, div_Ktau_yx_yy = divs
    grad_u2_x, grad_u2_y, grad_phi_x, grad_phi_y = grads

    # Equations
    momentum_x = div_phi_11_12 - div_delta_1_2_ast - S0
    momentum_y = div_phi_21_22 - T0
    kinetic_energy = div_phi_1_2_ast - div_delta_1_2_ast - (delta_1_line * grad_u2_x + delta_2_line * grad_u2_y) - 2 * C_D
    lateral_curvature = div_theta_1_2_o + (phi_1_ast * grad_phi_x + phi_2_ast * grad_phi_y) + 0.5 * (delta_1_line * grad_u2_y - delta_2_line * grad_u2_x) - (delta_1_o * grad_u2_x + delta_2_o * grad_u2_y) + C_D_x - 2 * C_D_o
    shear_stress_x = div_Ktau_xx_xy - Sx
    shear_stress_y = div_Ktau_yx_yy - Sy

    return [momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y]

def obj_func(equations: Iterable[float]) -> float:

    obj = 0.0

    for eq in equations:
        obj = obj + fabs(eq)

    return obj * obj

#----------------------------------------------#
# Tests
#----------------------------------------------#
def test_profiles():

    x = Parameters()
    freestream = Freestream()

    eta, U, W, dU_deta, dW_deta, R, S, T = profiles_func(x.to_list(), freestream.to_list())

    plt.figure()
    plt.plot(U, eta, label='U')
    plt.plot(W, eta, label='W')
    plt.legend()
    plt.grid()

    plt.figure()
    plt.plot(S, eta, label='S')
    plt.plot(T, eta, label='T')
    plt.legend()
    plt.grid()

    plt.figure()
    plt.plot(dU_deta, eta, label='dU_deta')
    plt.plot(dW_deta, eta, label='dW_deta')
    plt.legend()
    plt.grid()

    plt.figure()
    plt.plot(R, eta, label='R')
    plt.legend()
    plt.grid()

    plt.show()

    return

#----------------------------------------------#
# Gradient
#----------------------------------------------#
def find_minimum():

    x = Parameters()
    air = Freestream()
    divs = Divergentes()
    grads = Gradients()

    step = 1e-2

    for interaction in range(5):

        # Equations
        eq_ref = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())
        eq_delta = equations_func(x.convert_to_real(x.to_list_eps('delta')), divs.to_list(), grads.to_list(), air.to_list())
        eq_A = equations_func(x.convert_to_real(x.to_list_eps('A')), divs.to_list(), grads.to_list(), air.to_list())
        eq_B = equations_func(x.convert_to_real(x.to_list_eps('B')), divs.to_list(), grads.to_list(), air.to_list())
        eq_Psi = equations_func(x.convert_to_real(x.to_list_eps('Psi')), divs.to_list(), grads.to_list(), air.to_list())
        eq_Ctau1 = equations_func(x.convert_to_real(x.to_list_eps('Ctau1')), divs.to_list(), grads.to_list(), air.to_list())
        eq_Ctau2 = equations_func(x.convert_to_real(x.to_list_eps('Ctau2')), divs.to_list(), grads.to_list(), air.to_list())

        # Gradient
        grad = [
            - step * (obj_func(eq_delta) - obj_func(eq_ref)) / x.eps,
            - step * (obj_func(eq_A) - obj_func(eq_ref)) / x.eps,
            - step * (obj_func(eq_B) - obj_func(eq_ref)) / x.eps,
            - step * (obj_func(eq_Psi) - obj_func(eq_ref)) / x.eps,
            - step * (obj_func(eq_Ctau1) - obj_func(eq_ref)) / x.eps,
            - step * (obj_func(eq_Ctau2) - obj_func(eq_ref)) / x.eps,
        ]

        grad = [
            0,
            - step * (obj_func(eq_A) - obj_func(eq_ref)) / x.eps,
            0,
            0,
            0,
            0,
        ]

        # Increase
        x.add(grad)

        # New obj
        eq_ref2 = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())

        # Error
        print(grad)
        print(' {} - ({:e}; {:e}); [{:e}, {:e}, {:e}, {:e}, {:e}, {:e}]'.format(interaction, obj_func(eq_ref), obj_func(eq_ref2), eq_ref2[0], eq_ref2[1], eq_ref2[2], eq_ref2[3], eq_ref2[4], eq_ref2[5]))#x.delta, x.A, x.B, x.Psi, x.Ctau1, x.Ctau2))

    return

def plots():

    x = Parameters()
    air = Freestream()
    divs = Divergentes()
    grads = Gradients()

    # Delta
    delta = np.geomspace(1e-8, 1e-1, num=1000)
    obj_delta = np.zeros_like(delta)

    for i in range(obj_delta.size):
        x.delta = delta[i]
        eq = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())
        obj_delta[i] = obj_func(eq)
    
    # A
    A = np.linspace(-5, 5, num=1000)
    obj_A = np.zeros_like(A)

    for i in range(obj_A.size):
        x.A = A[i]
        eq = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())
        obj_A[i] = obj_func(eq)
    
    # B
    B = np.linspace(-5, 5, num=1000)
    obj_B = np.zeros_like(B)

    for i in range(obj_B.size):
        x.B = B[i]
        eq = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())
        obj_B[i] = obj_func(eq)
    
    # Psi
    Psi = np.linspace(-3.14, 3.14, num=1000)
    obj_Psi = np.zeros_like(Psi)

    for i in range(obj_Psi.size):
        x.Psi = Psi[i]
        eq = equations_func(x.convert_to_real(x.to_list()), divs.to_list(), grads.to_list(), air.to_list())
        obj_Psi[i] = obj_func(eq)

    plt.figure()
    plt.plot(delta, obj_delta)
    plt.grid()
    plt.title('delta')
    plt.xscale('log')

    plt.figure()
    plt.plot(A, obj_A)
    plt.grid()
    plt.title('A')
    
    plt.figure()
    plt.plot(B, obj_B)
    plt.grid()
    plt.title('B')

    plt.figure()
    plt.plot(Psi, obj_Psi)
    plt.grid()
    plt.title('Psi')

    plt.show()

    return

if __name__ == '__main__':
    plots()