import io
import os
import subprocess
import unittest
import yaml
import math
import numpy as np
import shutil
import uuid


def dihedral(p0, p1, p2, p3):
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


class Calculation:
    def __init__(self, input):
        self._input = input
        self._output = ''
        with open(self._input['systems'][0]['path']) as f:
            self._natoms = int(f.readline())
        self.energy = 0
        self.lowest_frequency = 0
        self.optimized_structure = []
        self.ts_guess_structure = []
        self.interpolated_trajectory = []
        self.optimized_trajectory = []
        self.tangent = []

    def run(self):
        self.__prepare_input()
        self.__execute()
        return self.__parse_output()

    def __prepare_input(self):
        with io.open('input.yaml', 'w', encoding='utf8') as outfile:
            yaml.dump(self._input, outfile, default_flow_style=False, allow_unicode=True)

    def __execute(self):
        # If the return code is non-zero, this automatically raises an exception
        self._output = subprocess.check_output(['readuct', 'input.yaml']).decode()

    def __parse_output(self):
        splitted_output = self._output.splitlines()
        success = False
        for index, line in enumerate(splitted_output):
            if "Completed all tasks, exiting!" in line:
                success = True
            if line.startswith('  The (electronic) energy is'):
                self.energy = float(splitted_output[index].split()[4])
            if line == '     #    cm^-1':
                self.lowest_frequency = float(splitted_output[index + 1].split()[1])
        if self._input['tasks'][0]['type'] in ['afir', 'geoopt', 'ts', 'irc']:
            self.optimized_structure = self.__load_xyz_file('default_system_output/default_system_output.xyz')
        if self._input['tasks'][0]['type'] in ['bspline']:
            self.interpolated_trajectory = self.__load_trajectory(
                'default_system_output/default_system_output_interpolated.xyz')
            if os.path.isfile('default_system_output/default_system_output_optimized.xyz'):
                self.optimized_trajectory = self.__load_trajectory(
                    'default_system_output/default_system_output_optimized.xyz')
            if os.path.isfile('default_system_output/default_system_output_tsguess.xyz'):
                self.ts_guess_structure = self.__load_xyz_file(
                    'default_system_output/default_system_output_tsguess.xyz')
            if os.path.isfile('default_system_output/default_system_output_tsguess-1.xyz'):
                self.ts_guess_structure_minus_1 = self.__load_xyz_file(
                    'default_system_output/default_system_output_tsguess-1.xyz')
            if os.path.isfile('default_system_output/default_system_output_tsguess+1.xyz'):
                self.ts_guess_structure_plus_1 = self.__load_xyz_file(
                    'default_system_output/default_system_output_tsguess+1.xyz')
            if os.path.isfile('default_system_output/tangent.dat'):
                f = open('default_system_output/tangent.dat', 'r')
                self.tangent = f.readlines()
                f.close()

        return success

    def __load_xyz_file(self, molecule_file):
        with open(molecule_file, "r") as f:
            coordinate_string = f.readlines()[2:]

        coordinates = []
        for coordinate in coordinate_string:
            split_coord = coordinate.split()
            element = split_coord[0]
            x = float(split_coord[1])
            y = float(split_coord[2])
            z = float(split_coord[3])
            coordinates.append([element, x, y, z])

        return coordinates

    def __load_trajectory(self, filename):
        lines_per_file = self._natoms + 2
        smallfile = None
        small_filename = str(uuid.uuid4())
        trajectory = []
        with open(filename) as trajectory_file:
            for lineno, line in enumerate(trajectory_file):
                if lineno % lines_per_file == 0:
                    if smallfile:
                        smallfile.close()
                        trajectory.append(self.__load_xyz_file(small_filename))
                    smallfile = open(small_filename, 'w')
                smallfile.write(line)
            if smallfile:
                smallfile.close()
                trajectory.append(self.__load_xyz_file(small_filename))
        os.remove(small_filename)

        return trajectory


class TestReaductFast(unittest.TestCase):
    def tearDown(self):
        if os.path.exists('default_system'):
            shutil.rmtree('default_system')
        if os.path.exists('default_system_output'):
            shutil.rmtree('default_system_output')
        if os.path.exists('default_system_output_2'):
            shutil.rmtree('default_system_output_2')
        if os.path.exists('default_system_ts'):
            shutil.rmtree('default_system_ts')
        if os.path.exists('input.yaml'):
            os.remove('input.yaml')

    def test_single_point_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6'
             }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'sp',
             }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        self.assertAlmostEqual(default_calculation.energy, -22.314858496)

        input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'AM1',
             'settings': {
                'molecular_charge': 1,
                'spin_multiplicity': 2,
                'unrestricted_calculation': 'true'
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'sp',
             }
        ]
        }

        calculation = Calculation(input)
        success = calculation.run()
        self.assertTrue(success)
        self.assertAlmostEqual(calculation.energy, -24.04220233)

    def test_hessian_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
                              'settings': {
                'self_consistence_criterion': 1e-8
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'hessian',
             }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        # The reference value was obtained form MOPAC as implemented in ADF 2016.107.
        self.assertAlmostEqual(default_calculation.lowest_frequency, -106.4, places=0)

        input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'DFTB2',
             'settings': {
                'molecular_charge': -1,
                'spin_multiplicity': 2,
                'unrestricted_calculation': 'true'
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'hessian',
             }
        ]
        }

        calculation = Calculation(input)
        success = calculation.run()
        self.assertTrue(success)
        self.assertAlmostEqual(calculation.lowest_frequency, -9073.1, places=0)

    def test_structure_optimization_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6'
             }
        ],
            'tasks': [
            {'input': ['default_system'],
             'output': ['default_system_output'],
             'type': 'geoopt',
             'settings': {'bfgs_use_gdiis': False}
             }
        ]
        }
        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = dihedral(np.array(loaded[2][1:]), np.array(loaded[0][1:]),
                          np.array(loaded[1][1:]), np.array(loaded[3][1:]))
        expected = 180.0
        # Due to the form of the PES in this case the fidelity of the
        #  test is only the first digit.
        self.assertAlmostEqual(abs(actual-expected), 0.0e0, places=1)

    def test_transition_state_optimization_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
                              'settings': {
                'self_consistence_criterion': 1e-8
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'output': ['default_system_output'],
             'type': 'ts',
             'settings': {'optimizer': 'ev'},
             },
            {'input': ['default_system_output'],
             'type': 'hessian'
             }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = dihedral(np.array(loaded[2][1:]), np.array(loaded[0][1:]),
                          np.array(loaded[1][1:]), np.array(loaded[3][1:]))
        expected = 0.0
        # Due to the form of the PES in this case the fidelity of the
        #  test is only the first digit.
        self.assertAlmostEqual(abs(actual-expected), 0.0e0, places=1)
        # The reference value was obtained form MOPAC as implemented in ADF 2016.107.
        self.assertAlmostEqual(default_calculation.lowest_frequency, -543.5, places=0)

    def test_transition_state_optimization_task_dimer(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {
                 'self_consistence_criterion': 1e-8
             }
             }
        ],
            'tasks': [
                {'input': ['default_system'],
                 'output': ['default_system_output'],
                 'type': 'ts',
                 'settings': {'optimizer': 'DIMER',
                              'geoopt_transform_coordinates': True,
                              'convergence_step_max_coefficient': 5e-5,
                              'convergence_step_rms': 1e-4,
                              'convergence_gradient_max_coefficient': 1e-5,
                              'convergence_gradient_rms': 5e-6,
                              'convergence_delta_value': 5e-8,
                              'convergence_max_iterations': 1000,
                              'convergence_requirement': 4,
                              'dimer_calculate_hessian_once': True,
                              'dimer_decrease_rotation_gradient_threshold': True,
                              'dimer_gradient_interpolation': True,
                              'dimer_radius': 0.02,
                              'dimer_phi_tolerance': 2e-3,
                              'dimer_rotation_gradient_first': 2e-7,
                              'dimer_rotation_gradient_other': 2e-4,
                              'dimer_lowered_rotation_gradient': 2e-3,
                              'dimer_grad_rmsd_threshold': 2e-3,
                              'dimer_trust_radius': 0.2,
                              'dimer_default_translation_step': 1.1,
                              'dimer_max_rotations_first_cycle': 200,
                              'dimer_max_rotations_other_cycle': 150,
                              'dimer_interval_of_rotations': 4,
                              'dimer_cycle_of_rotation_gradient_decrease': 6,
                              'dimer_max_backtracking': 6
                              },
                 },
                {'input': ['default_system_output'],
                 'type': 'hessian'
                 }
            ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = dihedral(np.array(loaded[2][1:]), np.array(loaded[0][1:]),
                          np.array(loaded[1][1:]), np.array(loaded[3][1:]))
        expected = 0.0
        # Due to the form of the PES in this case the fidelity of the
        #  test is only the first digit.
        self.assertAlmostEqual(abs(actual-expected), 0.0e0, places=1)
        # The reference value was obtained form MOPAC as implemented in ADF 2016.107.
        self.assertAlmostEqual(default_calculation.lowest_frequency, -543.8, places=1)

    def test_afir_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'br_2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
                              'settings': {
                'molecular_charge': -2
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'output': ['default_system_output'],
             'type': 'afir',
             'settings': {'allow_unconverged' : 'True'
                          },
             }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        expected_structure = [['Br', -6.0715944861, 0.0, 0.0], ['Br', 6.0715944861, 0.0, 0.0]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index])

    def test_irc_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2_ts.xyz'),
             'name': 'default_system_ts',
             'method_family': 'PM6',
                              'settings': {
                'self_consistence_criterion': 1e-8
            }
            }
        ],
            'tasks': [
            {'input': ['default_system_ts'],
             'output': ['default_system_output', 'default_system_output_2'],
             'type': 'irc',
             'settings': {'convergence_max_iterations': 150,
                          'convergence_delta_value': 1.0e-8,
                          'convergence_gradient_rms': 1.0e-7,
                          'convergence_step_rms': 1.0e-7,
                          'optimizer': 'lbfgs',
                          'lbfgs_linesearch': 'false',
                          'irc_initial_step_size': 1.0,
                          },
             },
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = dihedral(np.array(loaded[2][1:]), np.array(loaded[0][1:]),
                          np.array(loaded[1][1:]), np.array(loaded[3][1:]))
        expected = 180.0
        # Due to the form of the PES in this case the fidelity of the
        #  test is only the first digit.
        self.assertAlmostEqual(abs(actual-expected), 0.0e0, places=1)

    def test_bsplines_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'start.xyz'),
             'name': 'start',
             'method_family': 'PM6',
                              'settings': {
                'molecular_charge': -1
            }
            },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'end.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {
                'molecular_charge': -1
            }
            }
        ],
            'tasks': [
            {'input': ['start', 'end'],
             'output': ['default_system_output'],
             'type': 'bspline',
             'settings': {
                'extract_ts_guess': True,
                'convergence_max_iterations': 500,
                'tangent_file': "tangent.dat"
            },
            },
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        self.assertTrue(default_calculation.tangent)

        expected_interpolated_structure_5 = [['C', -2.1242184803, 1.9447120704, -0.7325866338],
                                             ['H', -1.3178297538, 2.3326493335, -0.1077572937],
                                             ['H', -3.0933405168, 2.1499835820, -0.2791638336],
                                             ['S', -2.0069401132, 2.6032604720, -2.4206019707],
                                             ['H', -1.9941622498, 0.8653072583, -0.7922433370],
                                             ['I', -3.3018395895, 9.0184385717, -2.7246106738],
                                             ['C', -2.5559097834, 5.3414727878, -2.6017902431],
                                             ['H', -0.4642288832, 5.8001953858, -2.4428873984],
                                             ['H', -3.8505950032, 5.6252075133, -0.8952195322],
                                             ['H', -2.9649581461, 4.4456368244, -4.4555935333],
                                             ['C', -3.8806282656, 5.3102305048, -1.9480106828],
                                             ['H', -4.5788393440, 5.9520250021, -2.4598131175],
                                             ['H', -4.2882251545, 4.3252001053, -1.9873623407],
                                             ['C', -2.4960576739, 5.3262409114, -4.0708497373],
                                             ['H', -3.0467094997, 6.1490768553, -4.4761105782],
                                             ['H', -1.4714523199, 5.3909419203, -4.4565570508],
                                             ['C', -1.3704237926, 5.7145090364, -1.8223057763],
                                             ['H', -1.5201445554, 6.6547134863, -1.3179918625],
                                             ['H', -1.1884174644, 4.9832220092, -1.0588847536]]
        for expected, actual in zip(expected_interpolated_structure_5, default_calculation.interpolated_trajectory[5]):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=3)

        expected_optimized_structure_0 = [['C', -2.0452993782, 1.5060589567, -0.6058933748],
                                          ['H', -1.2472703097, 1.9485506223, -0.0107946043],
                                          ['H', -3.0180758444, 1.7824248543, -0.2050791936],
                                          ['S', -1.8939176883, 2.1079356524, -2.3145874117],
                                          ['H', -1.9330206160, 0.4289364545, -0.5877846140],
                                          ['I', -3.0293834410, 8.3357711082, -2.9086214166],
                                          ['C', -2.6810543530, 6.0162055829, -2.6837688773],
                                          ['H', -0.5126814787, 6.0124812888, -2.4610926079],
                                          ['H', -3.9702666455, 5.9259625779, -0.9328618767],
                                          ['H', -2.9377848937, 4.5381453091, -4.2482732753],
                                          ['C', -3.9469126455, 5.6189712379, -1.9875730009],
                                          ['H', -4.8344619904, 6.0300161373, -2.4865910456],
                                          ['H', -4.0693736956, 4.5114056597, -1.9974566561],
                                          ['C', -2.5475866101, 5.5706429795, -4.1126930423],
                                          ['H', -3.1308458237, 6.2088090273, -4.7833723578],
                                          ['H', -1.5069871898, 5.5594830825, -4.4666781863],
                                          ['C', -1.4321618951, 5.9783403076, -1.8589462206],
                                          ['H', -1.3996069016, 6.8102467725, -1.1424972561],
                                          ['H', -1.3782291890, 5.0426360190, -1.2557753314]]
        for expected, actual in zip(expected_optimized_structure_0, default_calculation.optimized_trajectory[0]):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=3)

        expected_optimized_structure_5 = [['C', -2.1184451699, 1.9442058857, -0.7416076807],
                                          ['H', -1.3205787145, 2.3395799689, -0.1080666133],
                                          ['H', -3.0877375029, 2.1540241341, -0.2856069420],
                                          ['S', -2.0435142884, 2.7327329689, -2.4146150756],
                                          ['H', -1.9961763168, 0.8527042738, -0.8008708773],
                                          ['I', -3.2863516258, 8.9503428877, -2.7285828904],
                                          ['C', -2.5426181744, 5.2930436310, -2.5911348912],
                                          ['H', -0.4544312981, 5.8052050683, -2.4349633757],
                                          ['H', -3.8554347919, 5.6339385803, -0.8929320407],
                                          ['H', -2.9655710137, 4.4054545980, -4.4770763302],
                                          ['C', -3.8763702182, 5.3024729254, -1.9455029938],
                                          ['H', -4.5964840010, 5.9778374666, -2.4614119093],
                                          ['H', -4.3035406830, 4.2869443038, -1.9783021804],
                                          ['C', -2.4986221267, 5.3175885644, -4.0710860315],
                                          ['H', -3.0574857619, 6.1873612820, -4.4949775528],
                                          ['H', -1.4694381467, 5.3847494082, -4.4654922988],
                                          ['C', -1.3599614333, 5.7119312045, -1.8138109432],
                                          ['H', -1.5176291405, 6.6931607994, -1.3097624379],
                                          ['H', -1.1645301815, 4.9597456794,  -1.0345372847]]
        for expected, actual in zip(expected_optimized_structure_5, default_calculation.optimized_trajectory[5]):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_optimized_structure_9 = [['C', -2.1873537620, 2.2956345614, -0.8339412409],
                                          ['H', -1.3742773090, 2.6399283024, -0.1853274452],
                                          ['H', -3.1535522547, 2.4440305642, -0.3384315455],
                                          ['S', -2.0973580532, 2.9995203277, -2.5054136178],
                                          ['H', -2.0430755569, 1.2144039013, -0.9558103154],
                                          ['I', -3.5198045083, 9.5645725425, -2.5774020796],
                                          ['C', -2.4557941277, 4.8016865517, -2.5362073357],
                                          ['H', -0.4254668068, 5.6303666634, -2.4283232308],
                                          ['H', -3.7548576894, 5.3846034616, -0.8651056565],
                                          ['H', -2.9866967480, 4.3716300366, -4.6214497398],
                                          ['C', -3.8276007617, 5.0632379184, -1.9163608284],
                                          ['H', -4.3743412269, 5.8896320940, -2.4383907750],
                                          ['H', -4.4633063217, 4.1762356618, -1.9792868883],
                                          ['C', -2.4548345249, 5.1307192570, -4.0373750933],
                                          ['H', -2.9794004404, 6.1012911177, -4.2303011545],
                                          ['H', -1.4430244239, 5.2561089906, -4.4484601424],
                                          ['C', -1.3210333107, 5.5034440195, -1.7929934209],
                                          ['H', -1.6165746784, 6.5302868573, -1.4583875477],
                                          ['H', -1.0365680848, 4.9356908013, -0.9013722914]]
        for expected, actual in zip(expected_optimized_structure_9, default_calculation.optimized_trajectory[9]):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=3)

        expected_ts_guess_structure = [['C', -2.0916746662, 1.7935848435, -0.69528824374],
                                       ['H', -1.2959207114, 2.2089758649, -0.07356326423],
                                       ['H', -3.0636656621, 2.0287200143, -0.25860057205],
                                       ['S', -1.9990539050, 2.5353722331, -2.37828259150],
                                       ['H', -1.9754593759, 0.7038769724, -0.73035926397],
                                       ['I', -3.1926978350, 8.7144405623, -2.79027973720],
                                       ['C', -2.5915466959, 5.5577189166, -2.62222231602],
                                       ['H', -0.4689349351, 5.8755538451, -2.44251062809],
                                       ['H', -3.8962278558, 5.7348159996, -0.90390439036],
                                       ['H', -2.9534790989, 4.4348932937, -4.40042966651],
                                       ['C', -3.9011154561, 5.4144691850, -1.96020522079],
                                       ['H', -4.6903954056, 6.0017457704, -2.47101088188],
                                       ['H', -4.2202152208, 4.3471585970, -1.98310329607],
                                       ['C', -2.5178966222, 5.4083179011, -4.08639767871],
                                       ['H', -3.0888698983, 6.2064446325, -4.60724341309],
                                       ['H', -1.4804990712, 5.4419029240, -4.46929894629],
                                       ['C', -1.3820368438, 5.8077255163, -1.82684153651],
                                       ['H', -1.4719277940, 6.7453723392, -1.24273261215],
                                       ['H', -1.2333035361, 4.9719342194, -1.1080660913]]
        for expected, actual in zip(expected_ts_guess_structure, default_calculation.ts_guess_structure):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_tangent = [ -0.254487,
                             1.47209,
                             -0.456651,
                             -0.241359,
                             1.28398,
                             -0.332392,
                             -0.234254,
                             1.22696,
                             -0.265181,
                             -0.473553,
                             2.05341,
                             -0.342814,
                             -0.205065,
                             1.438,
                             -0.700538,
                             -0.893863,
                             2.21267,
                             0.600469,
                             0.480827,
                             -2.57208,
                             0.312053,
                             0.165121,
                             -0.686538,
                             0.078653,
                             0.393329,
                             -0.977861,
                             0.117168,
                             -0.112563,
                             -0.356626,
                             -0.766105,
                             0.250835,
                             -1.09098,
                             0.143841,
                             0.868093,
                             -0.198839,
                             0.0850303,
                             -0.813324,
                             -0.658956,
                             0.0558146,
                             0.185655,
                             -0.885504,
                             0.152927,
                             0.277957,
                             -0.128351,
                             1.04524,
                             0.115931,
                             -0.567032,
                             0.0210709,
                             0.226899,
                             -0.932896,
                             0.13167,
                             -0.432671,
                             -0.450127,
                             -0.619022,
                             0.696491,
                             -0.181327,
                             0.738766]

        for expected, actual in zip(expected_tangent, default_calculation.tangent):
                self.assertAlmostEqual(float(actual), expected, places=5)

    def test_bsplines_interpolation_task(self):
        input1 = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_1.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_2.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             }
        ],
            'tasks': [
            {'input': ['start', 'end'],
             'output': ['default_system_output'],
             'type': 'bspline',
             'settings': {
                'optimize': False,
                'extract_ts_guess': False
            },
            },
        ]
        }

        calculation1 = Calculation(input1)
        success = calculation1.run()
        self.assertTrue(success)

        self.assertEqual(len(calculation1.interpolated_trajectory), 10)

        # The entire structure is only translated in space; since this is removed during the interpolation, the
        # individual structures along the path remain the same.
        expected_interpolated_structure = [['H', 0.0, 0.0, 0.0], ['H', 1.0, 0.0, 0.0]]
        for structure in calculation1.interpolated_trajectory:
            for expected, actual in zip(expected_interpolated_structure, structure):
                for index in range(1, 4):
                    self.assertAlmostEqual(actual[index], expected[index], places=3)

        input2 = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_1.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_3.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             }
        ],
            'tasks': [
            {'input': ['start', 'end'],
             'output': ['default_system_output'],
             'type': 'bspline',
             'settings': {
                'optimize': False,
                'extract_ts_guess': False
            },
            },
        ]
        }

        calculation2 = Calculation(input2)
        success = calculation2.run()
        self.assertTrue(success)

        for index, structure in enumerate(calculation2.interpolated_trajectory):
            x = structure[0][1] - structure[1][1]
            y = structure[0][2] - structure[1][2]
            z = structure[0][3] - structure[1][3]
            bond_length = math.sqrt(x * x + y * y + z * z)
            self.assertAlmostEqual(bond_length, 1.0 + index * 1.0/9, places=3)
            self.assertAlmostEqual(structure[0][2], 0, places=3)
            self.assertAlmostEqual(structure[1][2], 0, places=3)
            self.assertAlmostEqual(structure[0][3], 0, places=3)
            self.assertAlmostEqual(structure[1][3], 0, places=3)


class TestReaductAll(TestReaductFast):
    def test_bsplines_task_2(self):
        input_2 = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'start2.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             'settings': {
                 'molecular_charge': -2
            }
            },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'end2.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {
                 'molecular_charge': -2
            }
            }
        ],
            'tasks': [
                {'input': ['start', 'end'],
                 'output': ['default_system_output'],
                 'type': 'bspline',
                 'settings': {
                     'extract_ts_guess': True,
                     'extract_ts_guess_neighbours': True
                },
                },
        ]
        }

        calculation_2 = Calculation(input_2)
        success = calculation_2.run()
        self.assertTrue(success)

        expected_interpolated_structure_5 = [['C', -0.2926301211, -0.5091949393, -0.5611081869],
                                             ['C', 0.3235912657, -1.2258151942, -1.6057373199],
                                             ['C', -0.0144474125, -2.5457401443, -1.8698819539],
                                             ['C', -0.9919141204, -3.1693804449, -1.0752366954],
                                             ['C', -1.6139086106, -2.4832194956, -0.0293627458],
                                             ['C', -1.2538310997, -1.1545867110, 0.2175253882],
                                             ['C', 0.0672617152, 0.8931254583, -0.3031703025],
                                             ['C', -0.9303353022, 1.8839025421, -0.2173441755],
                                             ['C', -0.6076179874, 3.2077170289, 0.0150302932],
                                             ['C', 0.7582311534, 3.5874176362, 0.1765690124],
                                             ['C', 1.7651483710, 2.5906575605, 0.0925846838],
                                             ['C', 1.4110263710, 1.2716629167, -0.1452425218],
                                             ['O', 0.9941919523, 4.8582018713, 0.3932096016],
                                             ['O', -1.2365094892, -4.4699293205, -1.4429642801],
                                             ['C', -2.2795874854, -5.2016442269, -0.7253635069],
                                             ['C', -2.3162792267, -6.5710107871, -1.4044317606],
                                             ['C', -3.5090029048, -7.3678151619, -0.9120834209],
                                             ['S', -3.4872416355, -9.0479750351, -1.6114472554],
                                             ['O', -2.1890593348, -9.6268848809, -1.2063583826],
                                             ['O', -3.6369327927, -8.8753328487, -3.0714187201],
                                             ['O', -4.6537558046, -9.6971309436, -0.9684141805],
                                             ['H', -3.5164456509, -7.4511750149, 0.1843190534],
                                             ['H', -4.4578970448, -6.8868889794, -1.1905159913],
                                             ['H', -1.3713122118, -7.1201536840, -1.2090526626],
                                             ['H', -2.3575524792, -6.4452106297, -2.5064881826],
                                             ['H', -3.2220282150, -4.6519261262, -0.8392612190],
                                             ['H', -1.9893397824, -5.2621434405, 0.3301629598],
                                             ['H', 0.4476401947, -3.1090278804, -2.6718513849],
                                             ['H', 1.0688609093, -0.7279007659, -2.2155019421],
                                             ['H', -1.7253466285, -0.6148784885, 1.0340737154],
                                             ['H', -2.3580037620, -2.9598570716, 0.5864966525],
                                             ['H', 2.1864838795, 0.5105097045, -0.2041005331],
                                             ['H', 2.8005820839, 2.8769509881, 0.2164245237],
                                             ['H', -1.3629986553, 3.9790049839, 0.0766561279],
                                             ['H', -1.9720151602, 1.5977133694, -0.3460339327],
                                             ['C', 4.2750594731, 5.9450478879, 1.7686352139],
                                             ['C', 2.7545155931, 5.8712827846, 1.8773755440],
                                             ['C', 2.2691403259, 6.7275508311, 3.0126267903],
                                             ['C', 2.0934979083, 6.2554404690, 0.5698727949],
                                             ['Br', 3.1679751728, 8.8878001902, 0.3844879305],
                                             ['H', 2.5917508777, 5.8534450730, -0.3037925720],
                                             ['H', 1.4855142823, 6.7760948564, 0.5286050390],
                                             ['H', 4.7629002115, 5.6550133631, 2.6769244502],
                                             ['H', 4.6338625713, 5.3135392951, 0.9860884211],
                                             ['H', 4.6090680335, 6.9530609510, 1.5499208894],
                                             ['H', 2.4453489380, 4.8202556586, 2.0644361707],
                                             ['H', 2.6902279949, 6.4144206361, 3.9512962793],
                                             ['H', 2.5597337878, 7.7595400312, 2.8674422839],
                                             ['H', 1.1843798515, 6.6854661284, 3.0754000101]]

        for expected, actual in zip(expected_interpolated_structure_5, calculation_2.interpolated_trajectory[5]):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=3)

        expected_optimized_structure_5 = [['C', -0.2915936065, -0.5094986883, -0.5623174051],
                                          ['C', 0.3264678238, -1.2256493765, -1.6082307988],
                                          ['C', -0.0114209354, -2.5455208297, -1.8724147897],
                                          ['C', -0.9928649173, -3.1715608806, -1.0753700819],
                                          ['C', -1.6168660954, -2.4844581659, -0.0271666415],
                                          ['C', -1.2566368011, -1.1557369860, 0.2199040081],
                                          ['C', 0.0708781156, 0.8942373580, -0.3032278977],
                                          ['C', -0.9294320833, 1.8872264604, -0.2167400712],
                                          ['C', -0.6009857912, 3.2106295370, 0.0164181388],
                                          ['C', 0.7724704580, 3.5980723514, 0.1760073293],
                                          ['C', 1.7697548847, 2.5983351494, 0.0946943413],
                                          ['C', 1.4114555101, 1.2720596050, -0.1445265864],
                                          ['O', 0.9734075716, 4.8794057807, 0.4016253252],
                                          ['O', -1.2351194897, -4.4708848092, -1.4449540686],
                                          ['C', -2.2810494882, -5.2027694781, -0.7251717901],
                                          ['C', -2.3165148871, -6.5721708560, -1.4050044633],
                                          ['C', -3.5098478333, -7.3678884805, -0.9119485523],
                                          ['S', -3.4878665536, -9.0485867356, -1.6116693744],
                                          ['O', -2.1860826368, -9.6283111029, -1.2051492489],
                                          ['O', -3.6372817330, -8.8745800007, -3.0749288376],
                                          ['O', -4.6570907314, -9.6983114335, -0.9669120284],
                                          ['H', -3.5171243557, -7.4516075524, 0.1848406711],
                                          ['H', -4.4590366039, -6.8869046762, -1.1908081199],
                                          ['H', -1.3701871386, -7.1216129668, -1.2087401529],
                                          ['H', -2.3579120487, -6.4456875203, -2.5085614553],
                                          ['H', -3.2262165485, -4.6504478085, -0.8400538952],
                                          ['H', -1.9896100637, -5.2630600713, 0.3339880332],
                                          ['H', 0.4504257228, -3.1105225826, -2.6747191771],
                                          ['H', 1.0736729010, -0.7264475865, -2.2192389027],
                                          ['H', -1.7285293206, -0.6144423871, 1.0368986199],
                                          ['H', -2.3634543436, -2.9623751267, 0.5902414834],
                                          ['H', 2.1879754825, 0.5094836776, -0.2037834976],
                                          ['H', 2.8039935231, 2.8770106314, 0.2168608773],
                                          ['H', -1.3625168654, 3.9799191651, 0.0768439301],
                                          ['H', -1.9727257453, 1.5977066512, -0.3462501159],
                                          ['C', 4.2881714749, 5.9408784068, 1.7710601147],
                                          ['C', 2.7746775073, 5.8589908996, 1.8876982457],
                                          ['C', 2.2610287309, 6.7328598666, 3.0226042217],
                                          ['C', 2.1462301371, 6.1456766578, 0.5625513308],
                                          ['Br', 3.1617517040, 8.8653857085, 0.3797540541],
                                          ['H', 2.6294460145, 5.8656436457, -0.3426105623],
                                          ['H', 1.3519026853, 6.8731851349, 0.5171604343],
                                          ['H', 4.7712593604, 5.6493946954, 2.6939143355],
                                          ['H', 4.6460722977, 5.3029484988, 0.9682164556],
                                          ['H', 4.6129982243, 6.9702604731, 1.5463983785],
                                          ['H', 2.4478026390, 4.8010804530, 2.0751746659],
                                          ['H', 2.6930172417, 6.4102422311, 3.9679678590],
                                          ['H', 2.5599671876, 7.7793243955, 2.8677869285],
                                          ['H', 1.1731394197, 6.6890786680, 3.0818887326]]

        for expected, actual in zip(expected_optimized_structure_5, calculation_2.optimized_trajectory[5]):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_structure = [['C', -0.2838885534, -0.5641543899, -0.5602777497],
                                       ['C', 0.3512566457, -1.2950378387, -1.5864473268],
                                       ['C', 0.0007975989, -2.6110318367, -1.8555936655],
                                       ['C', -1.0114412724, -3.2183977937, -1.0840727589],
                                       ['C', -1.6552090565, -2.5183904832, -0.0575306990],
                                       ['C', -1.2813530249, -1.1940438503, 0.1954402131],
                                       ['C', 0.0946406052, 0.8321622807, -0.2951511132],
                                       ['C', -0.8908210400, 1.8413220635, -0.2160782230],
                                       ['C', -0.5485531718, 3.1567876971, 0.0207038561],
                                       ['C', 0.8288354491, 3.5355403461, 0.1962467501],
                                       ['C', 1.8153023377, 2.5089276226, 0.1218207271],
                                       ['C', 1.4405663150, 1.1930401625, -0.1227268698],
                                       ['O', 1.0740580959, 4.7800235049, 0.4247469318],
                                       ['O', -1.2617290448, -4.5181410148, -1.4560074080],
                                       ['C', -2.3297844510, -5.2344637807, -0.7558245597],
                                       ['C', -2.3555529210, -6.6122829959, -1.4197262791],
                                       ['C', -3.5563734319, -7.4027391198, -0.9358740386],
                                       ['S', -3.5357524310, -9.0851225757, -1.6286939613],
                                       ['O', -2.2491938450, -9.6781685507, -1.1931716738],
                                       ['O', -3.6532805774, -8.9173389602, -3.0958894055],
                                       ['O', -4.7250135089, -9.7215406469, -1.0062698282],
                                       ['H', -3.5741082909, -7.4816625429, 0.1612715602],
                                       ['H', -4.5008865982, -6.9184776214, -1.2252607872],
                                       ['H', -1.4122439227, -7.1586386870, -1.2023963365],
                                       ['H', -2.3797759971, -6.4985858842, -2.5251922456],
                                       ['H', -3.2684967707, -4.6783460659, -0.9021601427],
                                       ['H', -2.0679959082, -5.2841892310, 0.3117795152],
                                       ['H', 0.4765822647, -3.1853219797, -2.6426139075],
                                       ['H', 1.1223901336, -0.8093780698, -2.1781558173],
                                       ['H', -1.7675904399, -0.6426461264, 0.9968964359],
                                       ['H', -2.4269267167, -2.9816277462, 0.5395917196],
                                       ['H', 2.2051612614, 0.4178518936, -0.1761590025],
                                       ['H', 2.8500479751, 2.7766911950, 0.2552085867],
                                       ['H', -1.2964524430, 3.9376221484, 0.0745622857],
                                       ['H', -1.9363870139, 1.5671097586, -0.3560876409],
                                       ['C', 4.3398165483, 6.1369121594, 1.7414086230],
                                       ['C', 2.8307504653, 6.0226573761, 1.8812529245],
                                       ['C', 2.3207059474, 6.8167341558, 3.0738519413],
                                       ['C', 2.1763108126, 6.3770880985, 0.5947797041],
                                       ['Br', 2.9274911236, 8.8731525371, 0.3880542131],
                                       ['H', 2.5831548449, 6.0462943988, -0.3335957360],
                                       ['H', 1.2233084516, 6.8620280426, 0.5987648730],
                                       ['H', 4.8447911921, 5.7882613941, 2.6325167942],
                                       ['H', 4.6883129553, 5.5551621111, 0.8905096543],
                                       ['H', 4.6464648032, 7.1802837985, 1.5775932374],
                                       ['H', 2.5263265565, 4.9418059278, 2.0140800951],
                                       ['H', 2.7783158494, 6.4481032286, 3.9907171600],
                                       ['H', 2.5883932535, 7.8759945790, 2.9760712715],
                                       ['H', 1.2350289456, 6.7381713115, 3.1530881035]]

        for expected, actual in zip(expected_ts_guess_structure, calculation_2.ts_guess_structure):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_minus_1_structure = [['C',  -0.2835054824, -0.5683845084, -0.5601238695],
                                               ['C',   0.3536965781, -1.2996255501, -1.5854365892],
                                               ['C',   0.0018588275, -2.6159666637, -1.8546638109],
                                               ['C',  -1.0125287763, -3.2215672461, -1.0848759384],
                                               ['C',  -1.6585943647, -2.5213462978, -0.0592771621],
                                               ['C',  -1.2832169454, -1.1966275750,  0.1938796236],
                                               ['C',   0.0971413124,  0.8280356615, -0.2945977938],
                                               ['C',  -0.8877117163,  1.8393552025, -0.2158039072],
                                               ['C',  -0.5431339744,  3.1533313894,  0.0211966155],
                                               ['C',   0.8344600248,  3.5318356812,  0.1971301303],
                                               ['C',   1.8192355287,  2.5042042708,  0.1240697030],
                                               ['C',   1.4429232062,  1.1885721260, -0.1208647199],
                                               ['O',   1.0775940565,  4.7790821444,  0.4285979096],
                                               ['O',  -1.2634590095, -4.5218338620, -1.4571213172],
                                               ['C',  -2.3335291516, -5.2368773311, -0.7578931926],
                                               ['C',  -2.3583096904, -6.6153286002, -1.4208993108],
                                               ['C',  -3.5598315354, -7.4052790630, -0.9375543701],
                                               ['S',  -3.5391795585, -9.0880245572, -1.6300035648],
                                               ['O',  -2.2532498857, -9.6818448842, -1.1921479684],
                                               ['O',  -3.6544529147, -8.9201891749, -3.0978515930],
                                               ['O',  -4.7303168077, -9.7233098917, -1.0088022278],
                                               ['H',  -3.5782573695, -7.4838621777,  0.1597378049],
                                               ['H',  -4.5040794311, -6.9206829231, -1.2277395714],
                                               ['H',  -1.4149597466, -7.1615701339, -1.2018942528],
                                               ['H',  -2.3813823642, -6.5024178681, -2.5268003863],
                                               ['H',  -3.2723244673, -4.6800195651, -0.9066856850],
                                               ['H',  -2.0735855691, -5.2858840719,  0.3109908033],
                                               ['H',   0.4791243402, -3.1908848376, -2.6410874044],
                                               ['H',   1.1269553334, -0.8148342983, -2.1761628720],
                                               ['H',  -1.7711017166, -0.6445345985,  0.9947933770],
                                               ['H',  -2.4326630483, -2.9835645678,  0.5369237270],
                                               ['H',   2.2067103794,  0.4113467824, -0.1741159208],
                                               ['H',   2.8543459662,  2.7696841663,  0.2580104342],
                                               ['H',  -1.2915629795,  3.9348370733,  0.0744603843],
                                               ['H',  -1.9339073009,  1.5650294072, -0.3568012428],
                                               ['C',   4.3487104293,  6.1494610345,  1.7404657193],
                                               ['C',   2.8401167510,  6.0303021214,  1.8854636504],
                                               ['C',   2.3236565541,  6.8238257916,  3.0791240772],
                                               ['C',   2.1709074381,  6.3840906957,  0.5938947586],
                                               ['Br',  2.9093701426,  8.8669507123,  0.3871916784],
                                               ['H',   2.5899497896,  6.0611722805, -0.3448315787],
                                               ['H',   1.1977368505,  6.8717211647,  0.6028764961],
                                               ['H',   4.8521978022,  5.7966048469,  2.6319663173],
                                               ['H',   4.6943144244,  5.5707144925,  0.8813589739],
                                               ['H',   4.6496710399,  7.1983045281,  1.5791717037],
                                               ['H',   2.5327596394,  4.9467362240,  2.0131062814],
                                               ['H',   2.7845827565,  6.4499508292,  3.9960356122],
                                               ['H',   2.5900701234,  7.8868283694,  2.9839224147],
                                               ['H',   1.2367545119,  6.7424832519,  3.1596680543]]

        for expected, actual in zip(expected_ts_guess_minus_1_structure, calculation_2.ts_guess_structure_minus_1):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_plus_1_structure = [['C', -0.2835993751, -0.5677091075, -0.5601507617],
                                              ['C',  0.3533937985, -1.2987668727, -1.5857092809],
                                              ['C',  0.0017119631, -2.6151561480, -1.8548753721],
                                              ['C', -1.0123004778, -3.2209912637, -1.0847685209],
                                              ['C', -1.6581241122, -2.5209281202, -0.0588985277],
                                              ['C', -1.2829149588, -1.1961549330,  0.1941854090],
                                              ['C',  0.0968524452,  0.8288053439, -0.2946976779],
                                              ['C', -0.8881880848,  1.8399271994, -0.2158112528],
                                              ['C', -0.5437733856,  3.1540015943,  0.0211457823],
                                              ['C',  0.8337837453,  3.5326249223,  0.1968785051],
                                              ['C',  1.8186773595,  2.5053213725,  0.1237371091],
                                              ['C',  1.4425630099,  1.1895492553, -0.1211335473],
                                              ['O',  1.0763223776,  4.7803392024,  0.4283220989],
                                              ['O', -1.2631277104, -4.5212501448, -1.4569870740],
                                              ['C', -2.3329282525, -5.2364868663, -0.7575138128],
                                              ['C', -2.3578271103, -6.6148339666, -1.4207179823],
                                              ['C', -3.5592570813, -7.4048478898, -0.9372582537],
                                              ['S', -3.5385881437, -9.0875732189, -1.6297931668],
                                              ['O', -2.2524650559, -9.6812301452, -1.1922942919],
                                              ['O', -3.6542554542, -8.9196592059, -3.0975972013],
                                              ['O', -4.7294812182, -9.7230241895, -1.0083133528],
                                              ['H', -3.5775534213, -7.4834909102,  0.1600299414],
                                              ['H', -4.5035632300, -6.9202923977, -1.2273137455],
                                              ['H', -1.4144380444, -7.1611140993, -1.2019723335],
                                              ['H', -2.3811124045, -6.5017641544, -2.5265974324],
                                              ['H', -3.2718069453, -4.6796724706, -0.9059185526],
                                              ['H', -2.0726163430, -5.2856239139,  0.3112705663],
                                              ['H',  0.4788043487, -3.1899615835, -2.6414881406],
                                              ['H',  1.1263589726, -0.8138067231, -2.1766758388],
                                              ['H', -1.7706227048, -0.6441852931,  0.9952917678],
                                              ['H', -2.4318849407, -2.9833297304,  0.5375549501],
                                              ['H',  2.2064996799,  0.4124787606, -0.1744571986],
                                              ['H',  2.8537804721,  2.7709247241,  0.2575365048],
                                              ['H', -1.2923798336,  3.9353616431,  0.0744888287],
                                              ['H', -1.9343578242,  1.5654078063, -0.3566798272],
                                              ['C',  4.3480861543,  6.1470305286,  1.7408354163],
                                              ['C',  2.8394509343,  6.0282634034,  1.8855553017],
                                              ['C',  2.3229067740,  6.8227946945,  3.0785032569],
                                              ['C',  2.1706263369,  6.3810810235,  0.5934862867],
                                              ['Br', 2.9122610908,  8.8668317025,  0.3870828791],
                                              ['H',  2.5905706342,  6.0589568142, -0.3449951435],
                                              ['H',  1.1991294809,  6.8719758592,  0.6018565488],
                                              ['H',  4.8512990288,  5.7948793620,  2.6327477853],
                                              ['H',  4.6938076932,  5.5675812454,  0.8822961234],
                                              ['H',  4.6492621394,  7.1957298375,  1.5787816007],
                                              ['H',  2.5317923830,  4.9449710998,  2.0138765857],
                                              ['H',  2.7835317932,  6.4494766521,  3.9957756977],
                                              ['H',  2.5897188789,  7.8856588150,  2.9825835461],
                                              ['H',  1.2359746185,  6.7418804867,  3.1587957975]]

        for expected, actual in zip(expected_ts_guess_plus_1_structure, calculation_2.ts_guess_structure_plus_1):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

if __name__ == '__main__':
    unittest.main()
