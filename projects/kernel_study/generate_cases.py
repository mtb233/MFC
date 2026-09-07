#!/usr/bin/env python3
import math
import json
from pathlib import Path


def calculate_eos_param(rhoFun, cFun, cpFun, pFun, tFun):
    gammaFun = cFun**2/(cpFun*tFun) + 1
    pi_infFun = rhoFun * cpFun * tFun * (gammaFun-1)/gammaFun - pFun
    return gammaFun, pi_infFun


def write_case_json(
  case_directory : str,
  cell_count : int,
  kernel_shape : int,
  kernel_deviation_mode : int,
  kernel_extent : int,
  epsilonb : float
):
  # Reference values for nondimensionalization 
  x_0 = 1.0e-3        # length - m
  rho_0 = 1.0e3      # density - kg/m3
  c_0 = 1475.0        # speed of sound - m/s
  p_0 = rho_0*c_0*c_0 # pressure - Pa
  T_0 = 298.0         # temperature - K

  # Host properties (water)
  p_atm = 101325.0                # Atmospheric pressure - Pa
  c_host = 1475.0                 # speed of sound - m/s
  rho_host = 1000.0               # density - kg/m3
  T_host = 298.0                  # temperature - K
  cp_host = 4180.0                # specific heat - J/kg*K
  mu_host = 1.0e-3
  [gamma_host, pi_inf_host]  = calculate_eos_param(rho_host, c_host, cp_host, p_atm, T_host)   # specific heat ratio and stiffness - Pa
  
  # Acoustic source properties
  p_amp = 2.0*p_atm # Amplitude of the acoustic source - Pa
  freq = 150e3      # Source frequency - Hz
  wlen = c_host/freq   # Wavelength - m

  # Lagrangian bubble's properties
  R_uni = 8314.0      # Universal gas constant - J/kmol/K
  MW_g = 28.0         # Molar weight of the gas - kg/kmol
  MW_v = 18.0         # Molar weight of the vapor - kg/kmol
  gamma_g = 1.4       # Specific heat ratio of the gas
  gamma_v = 1.333     # Specific heat ratio of the vapor
  p_v = 2350           # Vapor pressure of the host - Pa
  cp_g = 1.0e3         # Specific heat of the gas - J/kg/K
  cp_v = 2.1e3        # Specific heat of the vapor - J/kg/K
  k_g = 0.025         # Thermal conductivity of the gas - W/m/K
  k_v = 0.02          # Thermal conductivity of the vapor - W/m/K
  diffVapor = 2.5e-5  # Diffusivity coefficient of the vapor - m2/s
  sigma = 0.069   # Surface tension of the bubble - N/m

  # Domain and time set up
  x_b = -7.0e-3
  x_e = 7.0e-3
  y_b = -3.5e-3
  y_e = 3.5e-3
  z_b = -3.5e-3
  z_e = 3.5e-3
  stopTime  = 30.e-06 # stop time - sec
  saveTime  = 5.0e-6   # save time - sec

  with open(f"{case_directory}/case.json", "w") as case_json:
    json.dump({
      # Logistics ================================================
      'run_time_info'                : 'T',
      # ==========================================================

      # Computational Domain Parameters ==========================
      'x_domain%beg'                 : x_b/x_0,
      'x_domain%end'                 : x_e/x_0,
      'y_domain%beg'                 : y_b/x_0,
      'y_domain%end'                 : y_e/x_0,
      'z_domain%beg'                 : z_b/x_0,
      'z_domain%end'                 : z_e/x_0,
      'stretch_x'                    : 'F',
      'stretch_y'                    : 'F',
      'stretch_z'                    : 'F',
      'm'                            : 2*cell_count, 
      'n'                            : cell_count,
      'p'                            : cell_count,
      'cfl_adap_dt'                  : 'T',
      'cfl_target'                   : 0.4,
      'n_start'                      : 0,
      't_save'                       : saveTime*(c_0/x_0),
      't_stop'                       : stopTime*(c_0/x_0),
      # ==========================================================

      # Simulation Algorithm Parameters ==========================
      'model_eqns'                   : '5eq',
      'time_stepper'                 : 'rk3',
      'num_fluids'                   : 1,
      'num_patches'                  : 1,
      'viscous'                      : 'T',
      'mpp_lim'                      : 'F',
      'weno_order'                   : 7,
      'weno_eps'                     : 1.0E-16,
      'mapped_weno'                  :'T',
      'riemann_solver'               : 'hllc',
      'wave_speeds'                  : 'direct',
      'avg_state'                    : 'arithmetic',
      'bc_x%beg'                     :-6,
      'bc_x%end'                     :-6,
      'bc_y%beg'                     :-6,
      'bc_y%end'                     :-6,
      'bc_z%beg'                     :-6,
      'bc_z%end'                     :-6,
      # ==========================================================

      # Acoustic source ==========================================
      'acoustic_source'              : 'T',
      'num_source'                   : 1,
      'acoustic(1)%support'          : 3,     # Planar wave
      'acoustic(1)%pulse'            : 1,     # Sine wave
      'acoustic(1)%npulse'           : 1,     # number of pulses
      'acoustic(1)%mag'              : p_amp/p_0,  # Wave amplitude
      'acoustic(1)%wavelength'       : wlen/x_0,  # Wave length
      'acoustic(1)%length'           : 2.0*(z_e-z_b)/x_0,
      'acoustic(1)%height'           : 2.0*(y_e-y_b)/x_0,
      'acoustic(1)%loc(1)'           : -5.0e-03/x_0,
      'acoustic(1)%loc(2)'           : 0.0,
      'acoustic(1)%loc(3)'           : 0.0,
      'acoustic(1)%dir'              : 0.0,
      'acoustic(1)%delay'            : 0.0,
      # ==========================================================


      # Formatted Database Files Structure Parameters ============
      'format'                       : 'silo',
      'precision'                    : 'double',
      'prim_vars_wrt'                : 'T',
      'parallel_io'                  : 'T',
      'lag_db_wrt'                   : 'T',
      # ==========================================================

      # Patch 1: Water (left) ====================================
      'patch_icpp(1)%geometry'       : 9,
      'patch_icpp(1)%x_centroid'     : 0.0,
      'patch_icpp(1)%y_centroid'     : 0.0,
      'patch_icpp(1)%z_centroid'     : 0.0,
      'patch_icpp(1)%length_x'       : 2.0*(x_e-x_b)/x_0,
      'patch_icpp(1)%length_y'       : 2.0*(y_e-y_b)/x_0,
      'patch_icpp(1)%length_z'       : 2.0*(z_e-z_b)/x_0,
      'patch_icpp(1)%vel(1)'         : 0.0,
      'patch_icpp(1)%vel(2)'         : 0.0,
      'patch_icpp(1)%vel(3)'         : 0.0,
      'patch_icpp(1)%pres'           : p_atm/p_0,
      'patch_icpp(1)%alpha_rho(1)'   : rho_host/rho_0,
      'patch_icpp(1)%alpha(1)'       : 1.0,
      # ==========================================================

      # Lagrangian Bubbles ===========================
      'bubbles_lagrange'                 : 'T',
      'adap_dt'                          : 'T',       #Strang splitting
      "thermal"                          : 3,
      "polytropic"                       : "F",
      'bubble_model'                     : "keller_miksis",
      'lag_params%nBubs_glb'             : 1,
      'lag_params%solver_approach'       : 2,
      'lag_params%cluster_type'          : 2,
      'lag_params%pressure_corrector'    : 'T',
      'lag_params%kernel_shape'          : kernel_shape,
      'lag_params%kernel_deviation_mode' : kernel_deviation_mode,
      'lag_params%kernel_extent'         : kernel_extent,
      'lag_params%heatTransfer_model'    : 'T',
      'lag_params%massTransfer_model'    : 'T',
      'lag_params%epsilonb'              : epsilonb,
      'lag_params%valmaxvoid'            : 0.9,
      'lag_params%write_bubbles'         : 'T',
      'lag_params%write_bubbles_stats'   : 'T',
      # ==========================================================
      # Bubble parameters
      "bub_pp%R0ref": 1.0,
      "bub_pp%p0ref": 1.0,
      "bub_pp%rho0ref": 1.0,
      "bub_pp%T0ref": 1.0,
      "bub_pp%ss": sigma / (rho_0 * x_0 * c_0 * c_0),
      "bub_pp%pv": p_v / p_0,
      "bub_pp%vd": diffVapor / (x_0 * c_0),
      "bub_pp%mu_l": mu_host / (rho_0 * x_0 * c_0),
      "bub_pp%gam_v": gamma_v,
      "bub_pp%gam_g": gamma_g,
      "bub_pp%cp_v": cp_v * (T_0 / (c_0 * c_0)),
      "bub_pp%cp_g": cp_g * (T_0 / (c_0 * c_0)),
      "bub_pp%k_v": k_v * (T_0 / (x_0 * rho_0 * c_0 * c_0 * c_0)),
      "bub_pp%k_g": k_g * (T_0 / (x_0 * rho_0 * c_0 * c_0 * c_0)),
      "bub_pp%M_v": MW_v,
      "bub_pp%M_g": MW_g,
      "bub_pp%R_v": (R_uni / MW_v) * (T_0 / (c_0 * c_0)),
      "bub_pp%R_g": (R_uni / MW_g) * (T_0 / (c_0 * c_0)),

      # Fluids Physical Parameters ===============================
      # Host medium
      'fluid_pp(1)%gamma'            : 1.0/(gamma_host-1.0),
      'fluid_pp(1)%pi_inf'           : gamma_host*(pi_inf_host/p_0)/(gamma_host-1.0),
      'fluid_pp(1)%Re(1)'            : 1.0/(mu_host/(rho_0*c_0*x_0)),
    },
    case_json
  )


def generate_case(
  case_directory : str,
  cell_count : int,
  kernel_shape : int,
  kernel_deviation_mode : int,
  kernel_extent : int,
  epsilonb : float,
  bubble_position_x,
  bubble_position_y,
  bubble_position_z
) :
  x_0 = 1.0e-3        # length - m
  case_directory = Path(case_directory)
  case_directory.mkdir(parents=True, exist_ok=True)
  if not (case_directory / "input/").exists():
    (case_directory / "input").mkdir()
  with (case_directory / "input/lag_bubbles.dat").open("w") as bubbles:
    bubbles.write(f"   {(bubble_position_x/x_0):.6e}   {(bubble_position_y/x_0):.6e}    {(bubble_position_z/x_0):.6e}    0.000000E+00    0.000000E+00    0.000000E+00    0.500000E-01    0.000000E+00")
  write_case_json(
    case_directory,
    cell_count,
    kernel_shape,
    kernel_deviation_mode,
    kernel_extent,
    epsilonb
  )


# Generate cases
generate_case("projects/kernel_study/cases/df_large",    35,     1, 0, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/df_medium",   70,     1, 0, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/df_small",    140,    1, 0, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/df_tiny",     280,    1, 0, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)

generate_case("projects/kernel_study/cases/mc3_large",   35,     1, 1, 3, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc3_medium",  70,     1, 1, 3, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc3_small",   140,    1, 1, 3, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc3_tiny",    280,    1, 1, 3, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)

generate_case("projects/kernel_study/cases/mc0_large",   35,     1, 1, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc0_medium",  70,     1, 1, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc0_small",   140,    1, 1, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/mc0_tiny",    280,    1, 1, 0, 1.0,   1.0e-3, 1.0e-3, 1.0e-3)

generate_case("projects/kernel_study/cases/k1_large",    35,     1, 2, 0, 2*math.sqrt(2/math.pi),    1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k1_medium",   70,     1, 2, 0, 2*math.sqrt(2/math.pi),    1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k1_small",    140,    1, 2, 0, 2*math.sqrt(2/math.pi),    1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k1_tiny",     280,    1, 2, 0, 2*math.sqrt(2/math.pi),    1.0e-3, 1.0e-3, 1.0e-3)

generate_case("projects/kernel_study/cases/k2_large",    35,     1, 2, 0, 3.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k2_medium",   70,     1, 2, 0, 3.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k2_small",    140,    1, 2, 0, 3.0,   1.0e-3, 1.0e-3, 1.0e-3)
generate_case("projects/kernel_study/cases/k2_tiny",     280,    1, 2, 0, 3.0,   1.0e-3, 1.0e-3, 1.0e-3)
