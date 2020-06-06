# implementation of the grand-potential phase-field model based on M.Plapp PRE 84,031601(2011)+D.Cogswell PRE92, 011301(R) (2015)+Z.Hong ACS EL 2018, 3, 7, 1737–1743
# Mechanical suppression based on linear elasticity 
# w is the chemical potential, eta is the phase-field, pot is the electric potential
# For Zn metal battery
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 200
  xmax =200
  ny=200
  ymax=200
[]

[Variables]
  [./w]
  [../]
  [./eta]
  [../]
  [./pot]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]
[Functions]
[./ic_func_eta]
  type = ParsedFunction
  value = 0.5*(1.0-1.0*tanh((x-2)))
[../]
[./ic_func_c]
  type = ParsedFunction
  value = 0
[../]
  [./ic_func_pot]
  type = ParsedFunction
  value = -0.1*(1.0-tanh((x-2)*2))
[../]
[]
[ICs]
  [./eta]
    variable = eta
    type = FunctionIC
    function = ic_func_eta
  [../]
  [./w]
    variable = w
     type = FunctionIC
     function = ic_func_c
  [../]
  [./pot]
            variable = pot
    type = FunctionIC
    function = ic_func_pot
  [../]
[]
[BCs]
[./Periodic]
    [./all]
      auto_direction = 'y'
    [../]
 [../]
[./left_eta]
  type = DirichletBC
  variable = 'eta'
  boundary = 'left'
  value = 1
[../]
[./right_eta]
  type = DirichletBC
  variable = 'eta'
  boundary = 'right'
  value = 0
[../]
 [./left_w]
  type = DirichletBC
  variable = 'w'
  boundary = 'left'
  value = 0
[../]
[./right_w]
  type = DirichletBC
  variable = 'w'
  boundary = 'right'
  value = 0.0
[../]
  [./left_pot]
  type = DirichletBC
  variable = 'pot'
  boundary = 'left'
  value = -0.2
[../]
[./right_pot]
  type = DirichletBC
  variable = 'pot'
  boundary = 'right'
  value = 0.0
 [../]
  [./Pressure]
        [./load]
        #Applies the pressure
        boundary = right
        factor = 0.0 # MPa
        disp_x = disp_x
        disp_y = disp_y
      [../]
    [../]
        [./left_x]
      type = PresetBC
      variable = disp_x
      boundary = 'left'
      value = 0
    [../]
  []
[Kernels]
  [./w_dot]
    type = SusceptibilityTimeDerivative
    variable = w
    f_name = chi
    args = 'w' # in this case chi (the susceptibility) is simply a constant
  [../]
   [./Diffusion1]
    type = MatDiffusion
    variable=w
    D_name=D
   [../]
  [./Diffusion2]
    type = Diff
    variable = w
    cv=eta
    Q_name = zz
    QM_name = DN
    cp=pot
  [../]

  [./Noisew]
    type = LangevinNoise
    variable = w
    amplitude=0.0
  [../]
  [./coupled_etadot]
    type = CoupledSusceptibilityTimeDerivative
    variable = w
    v = eta
    f_name = ft
    args = 'eta'
  [../]
 [./elecN]
   type = Electronutrality
   variable = pot
   cp=eta
   cv =w
   Q_name = Le1    #before phi
   QM_name=0.
  [../]

  [./TensorMechanics]
    use_displaced_mesh = true
  [../]

 [./coupled_pos]
    type = CoupledSusceptibilityTimeDerivative
    variable = pot
    v = eta
    f_name = ft2
    args = 'eta'
  [../]
  [./BV]
    type = Kinetics
    variable = eta
    f_name = G
    cp=pot
    cv=eta
  [../]
  # Anisotropic surface energy
  [./anisoACinterface1]
    type = ACInterfaceKobayashi1
    variable = eta
    mob_name = L
  [../]
  [./anisoACinterface2]
    type = ACInterfaceKobayashi2
    variable = eta
    mob_name = L
  [../]
  [./AC_bulk]
    type = AllenCahn
    variable = eta
    f_name = FF
    mob_name = L
  [../]
 [./Noiseeta]
    type = LangevinNoise
    variable = eta
    amplitude=0.04
  [../]
  [./e_dot]
    type = TimeDerivative
    variable = eta
  [../]
[]
[AuxVariables]
  [./sigma11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
[AuxKernels]
  [./matl_sigma11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = sigma11_aux
  [../]
  [./matl_sigma12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = sigma12_aux
  [../]
  [./matl_sigma22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = sigma22_aux
  [../]
[]
[Materials]
[./constants]
  type = GenericConstantMaterial
#ul is the free energy density in the liquid phase
  prop_names  = 'kappa_op  M0  ft2   S1    S2   L    Ls    B   es    el  zz  A     ul    us    AA  dv     Cs  eigen'
  prop_values = '0.2   368 0.021  1000000 4.64  0.625  0.001  2.0  -13.8  4.0527  0 1.0   0.017225 13.8   77.38 0.93 10.92 0.4'
#  M0 M1 Normalized diffusion coefficient for liquid and solid B normalized constant for nF
#  M0=10^12*D    B=zF/(1000*R*T)   length 1 um, time 1s, energy is normalized by RT
[../]

  [./liquid_GrandPotential]
    type = DerivativeParsedMaterial
    function = 'ul-A*log(1+exp((w-el)/A))'
    args = 'w '
    f_name = f1
    material_property_names = 'A ul el'
    #outputs = exodus
  [../]
  [./solid_GrandPotential]
    type = DerivativeParsedMaterial
    function = 'us-A*log(1+exp((w-es)/A))'
    args = 'w'
    f_name = f2
    material_property_names = 'A us es'
   # outputs = exodus
  [../]
  [./switching_function]
    type = SwitchingFunctionMaterial
    eta ='eta'
    h_order = HIGH
  [../]
  [./barrier_function]
    type = BarrierFunctionMaterial
    eta = eta
  [../]
  [./total_GrandPotential]
    type = DerivativeTwoPhaseMaterial
    args = 'w'
    eta = eta
    fa_name = f1
    fb_name = f2
    derivative_order = 2
    W = 2.0
  [../]
  [./coupled_eta_function]
    type = DerivativeParsedMaterial
    function = '-(cs*(1+dv)-cl)*dh'  # in this code cs=-cs h=eta dh=1
    args = ' w eta'
    f_name = ft
    material_property_names = 'dh:=D[h,eta] h(eta) dv cs:=D[f2,w] cl:=D[f1,w]'
    derivative_order = 1
  [../]
  [./stiffness_a]
  type = ComputeElasticityTensor
  base_name = phasea
  # lambda, mu values
  C_ijkl = '0.5  0.25'           #a is for electrolyte b is for electrode
  fill_method = symmetric_isotropic
  # See RankFourTensor.h for details on fill methods
[../]
[./strain_a]
  type = ComputeFiniteStrain
  base_name = phasea
[../]
[./stress_a]
  type = ComputeFiniteStrainElasticStress
  base_name = phasea
[../]
[./stiffness_b]
  type = ComputeElasticityTensor
  base_name = phaseb
  # Stiffness tensor lambda, mu values
  # Note that the two phases could have different stiffnesses.
  # Try reducing the precipitate stiffness (to '1 1') rather than making it oversized
      C_ijkl = '38.8 38.8'
  fill_method = symmetric_isotropic
 # material_property_names = ' lam1 mu1'
[../]
[./strain_b]
  type = ComputeFiniteStrain
  base_name = phaseb
  eigenstrain_names = eigenstrain
[../]
  [./eigenstrain_b]
  type = ComputeEigenstrain
  base_name = phaseb
  eigen_base = '0.4 0.0 0.0'   # 0.4=1-porosity
  eigenstrain_name = eigenstrain
[../]
[./stress_b]
  type = ComputeFiniteStrainElasticStress
  base_name = phaseb
[../]
[./global_stress]
  type = TwoPhaseStressMaterial
  base_A = phasea
  base_B = phaseb
[../]
 [./PF]
  type = ParsedMaterial     
  function = '-sigma11_aux*eigen*1000/(Cs*3*8.314)'  #stress unit is GPa
  args = 'sigma11_aux'
  f_name = PF
  material_property_names = 'eigen Cs'
  outputs = other
[../]
  [./susceptibility]
      type = DerivativeParsedMaterial
      function = '-d2F1*(1-h)-d2F2*h*(dv+1)'
      args = 'w'
      f_name = chi
      derivative_order = 1
      #outputs = exodus
      material_property_names = 'h(eta) dv d2F1:=D[f1,w,w] d2F2:=D[f2,w,w]'
    [../]
   [./Diffusion coefficient]
    type = ParsedMaterial    
    function = 'cl*(h-1)*(1-h)*M0'  #cl is the negative value of concentration
    f_name = D
     args = 'eta w'
    material_property_names = 'M0 cl:=D[f1,w] h(eta)'
    outputs = other
  [../]
   [./Free]
    type = DerivativeParsedMaterial
    f_name = FF
    material_property_names = 'B'
    args='eta'
    function = 'B*eta*eta*(1-eta)*(1-eta)'
    derivative_order = 1
   #outputs = exodus
  [../]

  [./Convection coefficient]
    type = ParsedMaterial     
    function = 'AA*cl*(h-1)*(1-h)*M0'  #c is -c
    #function= '1'
    args = 'eta w'
    f_name = DN
    material_property_names = 'M0 AA cl:=D[f1,w] h(eta)'
  [../]
    [./Bultervolmer]
        type = DerivativeParsedMaterial
        function = 'Ls*(exp(pot*AA/2.)+58.556*cl*(1-h)*exp(-pot*AA/2.-PF))*dh'
        args = 'pot eta w'
        f_name = G
       derivative_order = 1
        material_property_names = 'Ls dh:=D[h,eta] h(eta) cl:=D[f1,w] AA PF'
        outputs = other
      [../]
 [./eta]
    type = ParsedMaterial
    f_name = etao
    args='eta'
    function = 'eta'
  # outputs = exodus
  [../]
  [./concentration]
    type = ParsedMaterial
    f_name = c
    args='eta w'
    material_property_names = 'h(eta) dFl:=D[f1,w]'
    function = '-dFl*(1-h)'
   outputs = other
  [../]
  [./material]
    type = InterfaceOrientationMaterial
    op = eta
    anisotropy_strength = 0.04
    mode_number = 6
    eps_bar = 0.4
  [../]
  [./Le1]
  type = ParsedMaterial
  f_name = Le1
  args = 'eta'
  material_property_names = 'S1 S2 h(eta)'
  function = 'S1*h+S2*(1-h)'
   #outputs = exodu
[../]

[]
[GlobalParams]
  enable_jit = false           # We are having some trouble with JIT, just forget about it
  displacements = 'disp_x disp_y'
[]
[Postprocessors]
  [./ETA]
    type = ElementIntegralMaterialProperty
    mat_prop = etao
    execute_on = 'initial timestep_end'
  [../]


[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      21                  preonly       lu           2'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type =Newton
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
    dt=0.004
    end_time = 400
[]

[Outputs]
  exodus = false
  csv = true
  execute_on = 'TIMESTEP_END'
  [./other]        # creates input_other.e
     type = Exodus
     interval = 10
  [../]
[]
