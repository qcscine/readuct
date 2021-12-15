__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

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
        if self._input['tasks'][0]['type'] in ['afir', 'geoopt', 'ts', 'irc', 'nt']:
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


class TestReaductBase(unittest.TestCase):
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


class TestReaductFast(TestReaductBase):
    def test_wrong_input(self):
        # Unrecognized keyword on top level
        wrong_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7}
             }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'sp',
             }
        ],
            'wrong_keyword': []
        }

        calculation = Calculation(wrong_input)
        with self.assertRaises(subprocess.CalledProcessError) as context:
            success = calculation.run()

        # Unrecognized keyword in system block
        wrong_system_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'wrong_keyword': 'something',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7}
             }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'sp',
             }
        ]
        }

        calculation = Calculation(wrong_system_input)
        with self.assertRaises(subprocess.CalledProcessError) as context:
            success = calculation.run()

        # Unrecognized keyword in task block
        wrong_task_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7}
             }
        ],
            'tasks': [
            {'input': ['default_system'],
             'type': 'sp',
             'wrong_keyword': 'something',
             }
        ]
        }

        calculation = Calculation(wrong_task_input)
        with self.assertRaises(subprocess.CalledProcessError) as context:
            success = calculation.run()

        # Impossible to achieve setting
        wrong_system_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'max_scf_iterations': 2}
             }
        ],
            'tasks': [
                {'input': ['default_system'],
                 'type': 'sp',
                 }
            ]
        }

        calculation = Calculation(wrong_system_input)
        with self.assertRaises(subprocess.CalledProcessError) as context:
            success = calculation.run()

        # Impossible to achieve setting but do not stop on error
        wrong_system_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'max_scf_iterations': 2}
             }
        ],
            'tasks': [
                {'input': ['default_system'],
                 'type': 'sp',
                 'settings': {'stop_on_error': False}
                 }
            ]
        }

        calculation = Calculation(wrong_system_input)
        success = calculation.run()
        assert success

    def test_single_point_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7}
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
                'spin_mode': 'unrestricted',
                'self_consistence_criterion': 1e-7,
                'density_rmsd_criterion' : 1e-7
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
                'self_consistence_criterion': 1e-8,
                'density_rmsd_criterion' : 1e-8
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
                'spin_mode': 'unrestricted'
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
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7}
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
        actual_dihedral = dihedral(np.array(loaded[2][1:]), np.array(loaded[0][1:]),
                                   np.array(loaded[1][1:]), np.array(loaded[3][1:]))
        expected_dihedral = 180.0
        difference = abs(actual_dihedral) - expected_dihedral

        # Due to the form of the PES in this case the fidelity of the test is only the first digit.
        self.assertAlmostEqual(difference, 0.0e0, places=1)

    def test_single_atom_opt_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'br.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7,
                          'density_rmsd_criterion' : 1e-7,
                          'molecular_charge': -1}
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
        self.assertTrue(all(abs(number) < 1e-12 for number in loaded[0][1:4]))

    def test_transition_state_optimization_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
             'name': 'default_system',
             'method_family': 'PM6',
                              'settings': {
                'self_consistence_criterion': 1e-8,
                'density_rmsd_criterion' : 1e-8
            }
            }
        ],
            'tasks': [
            {'input': ['default_system'],
             'output': ['default_system_output'],
             'type': 'ts',
             'settings': {'optimizer': 'ev',
                          'automatic_mode_selection': [2, 3]},
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
                 'self_consistence_criterion': 1e-8,
                 'density_rmsd_criterion' : 1e-8
            }
            }
        ],
            'tasks': [
                {'input': ['default_system'],
                 'output': ['default_system_output'],
                 'type': 'ts',
                 'settings': {'optimizer': 'DIMER',
                              'geoopt_coordinate_system': 'internal',
                              'convergence_step_max_coefficient': 1e-6,
                              'convergence_step_rms': 1e-4,
                              'convergence_gradient_max_coefficient': 1e-5,
                              'convergence_gradient_rms': 5e-6,
                              'convergence_delta_value': 1e-8,
                              'convergence_max_iterations': 1000,
                              'convergence_requirement': 4,
                              'dimer_calculate_hessian_once': True,
                              'dimer_decrease_rotation_gradient_threshold': True,
                              'dimer_gradient_interpolation': False,
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
                              'dimer_lbfgs_memory': 6
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
            {
                'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'br_2.xyz'),
                'name': 'default_system',
                'method_family': 'PM6',
                'settings': {'molecular_charge': -2,
                             'self_consistence_criterion': 1e-7,
                             'density_rmsd_criterion' : 1e-7}
            }
        ],
            'tasks': [
                {
                    'input': ['default_system'],
                    'output': ['default_system_output'],
                    'type': 'afir',
                    'settings': {'stop_on_error': False,
                                 'afir_coordinate_system': 'cartesianWithoutRotTrans'}
                }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        expected_struct = [['Br', -6.0715944861, 0.0, 0.0],
                           ['Br',  6.0715944861, 0.0, 0.0]]
        expected_dist = np.linalg.norm(np.array(expected_struct[0][1:]) - np.array(expected_struct[1][1:]))
        actual_struct = default_calculation.optimized_structure
        actual_dist = np.linalg.norm(np.array(actual_struct[0][1:]) - np.array(actual_struct[1][1:]))
        self.assertAlmostEqual(actual_dist, expected_dist, places=2)

    def test_afir_task_2(self):
        default_input = {'systems': [
            {
                'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sn2_start.xyz'),
                'name': 'default_system',
                'method_family': 'dftb3',
                'settings': {
                    'molecular_charge': -1,
                    'self_consistence_criterion': 1e-7,
                    'density_rmsd_criterion' : 1e-7
                }
            }
        ],
            'tasks': [{
                'input': ['default_system'],
                'output': ['default_system_output'],
                'type': 'afir',
                'settings': {
                    'afir_lhs_list': [0],
                    'afir_rhs_list': [1],
                    'convergence_max_iterations': 1000,
                },
            }]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = np.linalg.norm(np.array(loaded[0][1:]) - np.array(loaded[1][1:]))
        self.assertAlmostEqual(actual, 1.60, places=2)
        actual = np.linalg.norm(np.array(loaded[5][1:]) - np.array(loaded[1][1:]))
        self.assertAlmostEqual(actual, 3.62, places=2)

    def test_afir_task_3(self):
        default_input = {'systems': [
            {
                'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'br_2_alternative.xyz'),
                'name': 'default_system',
                'method_family': 'test',
                'settings': {'self_consistence_criterion': 1e-7}
            }
        ],
            'tasks': [
                {
                    'input': ['default_system'],
                    'output': ['default_system_output'],
                    'type': 'afir',
                    'settings': {'afir_coordinate_system': 'cartesian'}
                }
        ]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        expected_structure = [['Br', -0.1437805410, 0.0, 0.0], ['Br', 2.1437805410, 0.0, 0.0]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=6)

    def test_nt_task(self):
        default_input = {'systems': [{
            'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sn2_start.xyz'),
            'name': 'default_system',
            'method_family': 'dftb3',
            'settings': {
                'molecular_charge': -1,
                'self_consistence_criterion': 1e-7,
                'density_rmsd_criterion' : 1e-7
            }
        }],
            'tasks': [{
                'input': ['default_system'],
                'output': ['default_system_output'],
                'type': 'nt',
                'settings': {
                    'nt_lhs_list': [0],
                    'nt_rhs_list': [1],
                    'nt_total_force_norm': 0.01,
                    'nt_movable_side': 'rhs',
                    'nt_constrained_atoms': [0],
                    'nt_coordinate_system': 'cartesian',
                    'sd_factor': 0.5,
                    'convergence_max_iterations': 1000,
                    'convergence_attractive_stop': 0.9,
                },
            }]
        }

        default_calculation = Calculation(default_input)
        success = default_calculation.run()
        self.assertTrue(success)
        loaded = default_calculation.optimized_structure
        actual = np.linalg.norm(np.array(loaded[0][1:]) - np.array(loaded[1][1:]))
        self.assertAlmostEqual(round(actual, 2), 2.40, places=2)
        actual = np.linalg.norm(np.array(loaded[5][1:]) - np.array(loaded[1][1:]))
        # The following check yields 2.3 in release mode, but 2.5 in debug mode
        self.assertAlmostEqual(actual, 2.4, places=0)

    def test_irc_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2_ts.xyz'),
             'name': 'default_system_ts',
             'method_family': 'PM6',
                              'settings': {
                'self_consistence_criterion': 1e-8,
                'density_rmsd_criterion' : 1e-8
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
                          'lbfgs_linesearch': False,
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
        self.assertAlmostEqual(abs(actual)-expected, 0.0e0, places=1)

    def test_bsplines_task(self):
        default_input = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'start.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             'settings': {
                'molecular_charge': -1,
                'self_consistence_criterion': 1e-7,
                'density_rmsd_criterion' : 1e-7
            }
            },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'end.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {
                'molecular_charge': -1,
                'self_consistence_criterion': 1e-7,
                'density_rmsd_criterion' : 1e-7
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

        expected_optimized_structure_5 = [['C', -2.1181172732, 1.9836714746, -0.7785151073],
                                          ['H', -1.3197130037, 2.3197032410, -0.1209773277],
                                          ['H', -3.0749952310, 2.1284188534, -0.2889455979],
                                          ['S', -2.1290732494, 2.9938291531, -2.3413334667],
                                          ['H', -1.9925978376, 0.8922940213, -0.8965371237],
                                          ['I', -3.2139748891, 8.6749700156, -2.7566662543],
                                          ['C', -2.5448866471, 5.2462475016, -2.5688437573],
                                          ['H', -0.4569520827, 5.8340346354, -2.3917776711],
                                          ['H', -3.8760443571, 5.6605399299, -0.9159226700],
                                          ['H', -2.9469637601, 4.3866439507, -4.5024463007],
                                          ['C', -3.8865163851, 5.2895380992, -1.9426096937],
                                          ['H', -4.5867750036, 5.9879381742, -2.4454349478],
                                          ['H', -4.3288978721, 4.2786954278, -1.9496806504],
                                          ['C', -2.5056090597, 5.2940482309, -4.0577615825],
                                          ['H', -3.0403495176, 6.1836946320, -4.4833171205],
                                          ['H', -1.5023268880, 5.3793837211, -4.4799397886],
                                          ['C', -1.3667703121, 5.7104540747, -1.7968727269],
                                          ['H', -1.5112735075, 6.7101927824, -1.3232094348],
                                          ['H', -1.1130837127, 4.9787257115, -1.0095491274]]
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

        expected_ts_guess_structure = [['C', -2.0893811632, 1.8112813643, -0.7196906764],
                                       ['H', -1.2935641112, 2.1829003261, -0.0805074231],
                                       ['H', -3.0517651721, 1.9975594487, -0.2569168881],
                                       ['S', -2.0752328593, 2.7583694233, -2.3038229825],
                                       ['H', -1.9706067739, 0.7259654079, -0.8094247706],
                                       ['I', -3.1182629953, 8.4416715568, -2.8184112541],
                                       ['C', -2.6015580940, 5.5631014596, -2.6144149818],
                                       ['H', -0.4724691439, 5.9029663299, -2.4013512827],
                                       ['H', -3.9194874337, 5.7603494829, -0.9246580995],
                                       ['H', -2.9281905810, 4.4248542374, -4.4059143804],
                                       ['C', -3.9130661109, 5.4109761915, -1.9597971721],
                                       ['H', -4.6994444017, 6.0068074456, -2.4549206314],
                                       ['H', -4.2184621408, 4.3467867027, -1.9607180093],
                                       ['C', -2.5278053613, 5.3997859411, -4.0799929203],
                                       ['H', -3.0779558916, 6.2020757469, -4.6168466792],
                                       ['H', -1.5141037367, 5.4402061390, -4.4829731198],
                                       ['C', -1.3906666201, 5.8138762167, -1.8132369058],
                                       ['H', -1.4531850006, 6.7608993053, -1.2392958378],
                                       ['H', -1.1997129981, 4.9825909046, -1.1074463345]]
        for expected, actual in zip(expected_ts_guess_structure, default_calculation.ts_guess_structure):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_tangent = [-0.239278,
                            1.55754,
                            -0.602502,
                            -0.238813,
                            1.21556,
                            -0.3764,
                            -0.177045,
                            1.14689,
                            -0.305006,
                            -0.710881,
                            2.84015,
                            -0.168863,
                            -0.196327,
                            1.49092,
                            -0.906734,
                            -0.721461,
                            1.53762,
                            0.510299,
                            0.528091,
                            -3.00172,
                            0.439134,
                            0.191196,
                            -0.570724,
                            0.158384,
                            0.345808,
                            -0.838374,
                            0.0702003,
                            -0.133971,
                            -0.49866,
                            -0.943062,
                            0.237056,
                            -1.12668,
                            0.158944,
                            0.923427,
                            -0.0650731,
                            0.0936308,
                            -1.02015,
                            -0.753986,
                            0.154201,
                            0.181972,
                            -0.988692,
                            0.195317,
                            0.306691,
                            -0.0226949,
                            1.11677,
                            0.0742686,
                            -0.569341,
                            -0.0316138,
                            0.240139,
                            -0.943181,
                            0.187367,
                            -0.489968,
                            -0.283512,
                            -0.708076,
                            0.899243,
                            -0.126046,
                            0.958014]
        for expected, actual in zip(expected_tangent, default_calculation.tangent):
            self.assertAlmostEqual(float(actual), expected, places=3)

    def test_bsplines_interpolation_task(self):
        input1 = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_1.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7, 'density_rmsd_criterion' : 1e-7}
             },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_2.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7, 'density_rmsd_criterion' : 1e-7}
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
             'settings': {'self_consistence_criterion': 1e-7, 'density_rmsd_criterion' : 1e-7}
             },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2_3.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {'self_consistence_criterion': 1e-7, 'density_rmsd_criterion' : 1e-7}
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


class TestReaductSlow(TestReaductBase):
    def test_bsplines_task_2(self):
        input_2 = {'systems': [
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'start2.xyz'),
             'name': 'start',
             'method_family': 'PM6',
             'settings': {
                 'molecular_charge': -2,
                 'self_consistence_criterion': 1e-7,
                 'density_rmsd_criterion' : 1e-7
            }
            },
            {'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'end2.xyz'),
             'name': 'end',
             'method_family': 'PM6',
             'settings': {
                 'molecular_charge': -2,
                 'self_consistence_criterion': 1e-7,
                 'density_rmsd_criterion' : 1e-7
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

        expected_optimized_structure_5 = [['C', -0.2913361885, -0.5085154474, -0.5623099448],
                                          ['C', 0.3252158514, -1.2257408421, -1.6073473033],
                                          ['C', -0.0118367883, -2.5449574878, -1.8721453307],
                                          ['C', -0.9931931891, -3.1718166238, -1.0750920962],
                                          ['C', -1.6156438366, -2.4833682044, -0.027618957],
                                          ['C', -1.2563085416, -1.1555762908, 0.2198290907],
                                          ['C', 0.0693719082, 0.893694128, -0.3032358827],
                                          ['C', -0.9306026038, 1.8850574493, -0.2171788829],
                                          ['C', -0.6055149806, 3.2094448925, 0.0158812959],
                                          ['C', 0.7693759786, 3.5963824308, 0.1765486546],
                                          ['C', 1.7676936569, 2.5956508761, 0.0939786932],
                                          ['C', 1.4107661486, 1.2707241932, -0.1451387537],
                                          ['O', 0.9795284802, 4.8655061807, 0.3968347178],
                                          ['O', -1.2351960708, -4.470005023, -1.4443620639],
                                          ['C', -2.2802692615, -5.2024038961, -0.7253213327],
                                          ['C', -2.3164004555, -6.5716345535, -1.4046627865],
                                          ['C', -3.5093461232, -7.3675560826, -0.9118639442],
                                          ['S', -3.4876243734, -9.0474673051, -1.6112776633],
                                          ['O', -2.1862750593, -9.6280277421, -1.2054085992],
                                          ['O', -3.6371823152, -8.8747374206, -3.0743038321],
                                          ['O', -4.656168907, -9.6981578259, -0.9670919422],
                                          ['H', -3.5166611017, -7.4513451563, 0.1846511208],
                                          ['H', -4.4583991615, -6.8868360059, -1.1906118477],
                                          ['H', -1.3704864294, -7.1208796036, -1.2088587615],
                                          ['H', -2.3576969927, -6.4452902712, -2.5077021708],
                                          ['H', -3.2245516533, -4.6508527379, -0.8395697308],
                                          ['H', -1.989180945, -5.2626020184, 0.3327041899],
                                          ['H', 0.4489805756, -3.1096562502, -2.6733807174],
                                          ['H', 1.0713114983, -0.7269223185, -2.2175793044],
                                          ['H', -1.7269081588, -0.614572543, 1.0356411733],
                                          ['H', -2.3608105418, -2.9611905127, 0.5886581366],
                                          ['H', 2.1873006944, 0.5098521332, -0.20406689],
                                          ['H', 2.8014246433, 2.8771414539, 0.2164770877],
                                          ['H', -1.3634133417, 3.9797229718, 0.0767055985],
                                          ['H', -1.9728261926, 1.5975641546, -0.3461968412],
                                          ['C', 4.2762843835, 5.9423740366, 1.7692146273],
                                          ['C', 2.767098269, 5.867506011, 1.8789608711],
                                          ['C', 2.2620980846, 6.7302541376, 3.0198654072],
                                          ['C', 2.1657947098, 6.1672566018, 0.5702394251],
                                          ['Br', 3.1649614199, 8.8760394196, 0.3823374628],
                                          ['H', 2.6091065425, 5.861068202, -0.3216152023],
                                          ['H', 1.3820763722, 6.8533610385, 0.5203619477],
                                          ['H', 4.7670859809, 5.6514052831, 2.6872665667],
                                          ['H', 4.640308791, 5.3060976603, 0.9750685209],
                                          ['H', 4.6115712273, 6.9638143597, 1.5475774429],
                                          ['H', 2.4452594773, 4.8096227003, 2.0681793154],
                                          ['H', 2.692473972, 6.4114628379, 3.9611313748],
                                          ['H', 2.5606727153, 7.7718622016, 2.8670770911],
                                          ['H', 1.1780718327, 6.6872468091, 3.0787509696]]
        for expected, actual in zip(expected_optimized_structure_5, calculation_2.optimized_trajectory[5]):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_structure = [['C', -0.2844856556, -0.5573491651, -0.5605121762],
                                       ['C', 0.3473335304, -1.2876032261, -1.5879528624],
                                       ['C', -0.0009076526, -2.6033390816, -1.8571697552],
                                       ['C', -1.0097894086, -3.2136879211, -1.0828845115],
                                       ['C', -1.6498442203, -2.5136268406, -0.0546729532],
                                       ['C', -1.2783759254, -1.1897570797, 0.1980401333],
                                       ['C', 0.090775981, 0.8384132721, -0.2960149614],
                                       ['C', -0.8959748215, 1.8446322737, -0.2165109201],
                                       ['C', -0.5575792063, 3.1619542994, 0.0198863326],
                                       ['C', 0.8210255333, 3.5411516082, 0.1948479243],
                                       ['C', 1.808674994, 2.5162541862, 0.1182526522],
                                       ['C', 1.4365913664, 1.1999903886, -0.1257286389],
                                       ['O', 1.0681742184, 4.7802941227, 0.4181404767],
                                       ['O', -1.2588981738, -4.512047223, -1.4542314582],
                                       ['C', -2.323788112, -5.2307184604, -0.752624991],
                                       ['C', -2.3512135422, -6.6074114401, -1.4178004767],
                                       ['C', -3.5508380213, -7.3985970462, -0.9331867829],
                                       ['S', -3.5303517034, -9.0801493485, -1.6264910486],
                                       ['O', -2.2425339723, -9.6724698612, -1.194723501],
                                       ['O', -3.6514402045, -8.9128362866, -3.0930156022],
                                       ['O', -4.7167733421, -9.7188593515, -1.0021726398],
                                       ['H', -3.5674848687, -7.4781722237, 0.1636404416],
                                       ['H', -4.4957365101, -6.9150053939, -1.2213512757],
                                       ['H', -1.4079806275, -7.1539201482, -1.2031957593],
                                       ['H', -2.377200989, -6.4924767473, -2.5225548297],
                                       ['H', -3.2623135951, -4.6757255561, -0.8949620316],
                                       ['H', -2.0590740232, -5.2814608241, 0.3129388391],
                                       ['H', 0.4723394246, -3.1763676842, -2.644804135],
                                       ['H', 1.1147987638, -0.8008544408, -2.1810006047],
                                       ['H', -1.7617728295, -0.6397216076, 1.0000271526],
                                       ['H', -2.4174634767, -2.9783883914, 0.5435564984],
                                       ['H', 2.2026893628, 0.4281143161, -0.1794301424],
                                       ['H', 2.8427466309, 2.7877336611, 0.2506942682],
                                       ['H', -1.3044463244, 3.9421476804, 0.0746937046],
                                       ['H', -1.9403810317, 1.5703776982, -0.3549657643],
                                       ['C', 4.3221532107, 6.1170461744, 1.7427851373],
                                       ['C', 2.8166160695, 6.0120498567, 1.8733146799],
                                       ['C', 2.3151761549, 6.8051036724, 3.0652760117],
                                       ['C', 2.1876033071, 6.3632910516, 0.5982506144],
                                       ['Br', 2.9566635981, 8.8870690252, 0.3895982106],
                                       ['H', 2.5697944221, 6.0224303283, -0.3141057114],
                                       ['H', 1.2679187381, 6.8451892144, 0.5926220494],
                                       ['H', 4.8326290609, 5.7753926047, 2.6327048926],
                                       ['H', 4.6781300042, 5.5308609907, 0.9056717687],
                                       ['H', 4.6412712411, 7.1513985713, 1.5751485943],
                                       ['H', 2.5157186723, 4.9349097757, 2.0148051292],
                                       ['H', 2.7685610836, 6.4453091633, 3.9815884408],
                                       ['H', 2.5860723729, 7.8582326785, 2.9633409738],
                                       ['H', 1.233190497, 6.7311987353, 3.1422386072]]
        for expected, actual in zip(expected_ts_guess_structure, calculation_2.ts_guess_structure):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_minus_1_structure = [['C', -0.2844479495, -0.5576205499, -0.5605019027],
                                               ['C', 0.34745597310, -1.2879468087, -1.5878447900],
                                               ['C', -0.0008475936, -2.6036634450, -1.8570860420],
                                               ['C', -1.0098812598, -3.2139199388, -1.0829278516],
                                               ['C', -1.6500337284, -2.5137948215, -0.0548235999],
                                               ['C', -1.2784978735, -1.1899466846, 0.19791857780],
                                               ['C', 0.09089432890, 0.83810615400, -0.2959748538],
                                               ['C', -0.8957823519, 1.84440748950, -0.2165072548],
                                               ['C', -0.5573131181, 3.16169067830, 0.01990846820],
                                               ['C', 0.82131106920, 3.54084165070, 0.19494962300],
                                               ['C', 1.80890193790, 2.51581220150, 0.11838713190],
                                               ['C', 1.43673488270, 1.19959775200, -0.1256208761],
                                               ['O', 1.06866917680, 4.77982261440, 0.41825853580],
                                               ['O', -1.2590301106, -4.5122806262, -1.4542859428],
                                               ['C', -2.3240296525, -5.2308755408, -0.7527766155],
                                               ['C', -2.3514068503, -6.6076099788, -1.4178733841],
                                               ['C', -3.5510683621, -7.3987694753, -0.9333052422],
                                               ['S', -3.5305889078, -9.0803310291, -1.6265755882],
                                               ['O', -2.2428470676, -9.6727163768, -1.1946644085],
                                               ['O', -3.6515193306, -8.9130479975, -3.0931188221],
                                               ['O', -4.7171093183, -9.7189740633, -1.0023677695],
                                               ['H', -3.5677670689, -7.4783211777, 0.16352368720],
                                               ['H', -4.4959437500, -6.9151618521, -1.2215219643],
                                               ['H', -1.4081890452, -7.1541034662, -1.2031643542],
                                               ['H', -2.3773092757, -6.4927387788, -2.5226370276],
                                               ['H', -3.2625227171, -4.6758639461, -0.8952695681],
                                               ['H', -2.0594622048, -5.2815654520, 0.31282847210],
                                               ['H', 0.47246883390, -3.1767380243, -2.6446450901],
                                               ['H', 1.11503968950, -0.8012652465, -2.1807969997],
                                               ['H', -1.7619660853, -0.6398613374, 0.99982901360],
                                               ['H', -2.4177774340, -2.9784835898, 0.54330553900],
                                               ['H', 2.20277465450, 0.42766053040, -0.1792933322],
                                               ['H', 2.84297597480, 2.78723716470, 0.25088430400],
                                               ['H', -1.3041187926, 3.94193887740, 0.07468253030],
                                               ['H', -1.9402006571, 1.57022677250, -0.3550144230],
                                               ['C', 4.32240771090, 6.11801680390, 1.74263827780],
                                               ['C', 2.81688790200, 6.01285260830, 1.87328348100],
                                               ['C', 2.31547257370, 6.80551877460, 3.06552640160],
                                               ['C', 2.18769936590, 6.36439590560, 0.59840547100],
                                               ['Br', 2.9555078120, 8.88713369180, 0.38963912180],
                                               ['H', 2.56957291320, 6.02332456270, -0.3140596609],
                                               ['H', 1.26731511690, 6.84512714440, 0.59302446120],
                                               ['H', 4.83299200100, 5.77608217120, 2.63239945630],
                                               ['H', 4.67833851590, 5.53211090400, 0.90528911160],
                                               ['H', 4.64143555070, 7.15243779660, 1.57530209820],
                                               ['H', 2.51611006690, 4.93560816140, 2.01450804390],
                                               ['H', 2.76898306610, 6.44549792370, 3.98169973630],
                                               ['H', 2.58621326670, 7.85870940440, 2.96387548450],
                                               ['H', 1.23349812240, 6.73144246960, 3.14259033590]]
        for expected, actual in zip(expected_ts_guess_minus_1_structure, calculation_2.ts_guess_structure_minus_1):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)

        expected_ts_guess_plus_1_structure = [['C', -0.2845233666, -0.5570777854, -0.5605224457],
                                              ['C',  0.3472110826, -1.2872596457, -1.5880609288],
                                              ['C', -0.0009677203, -2.6030147215, -1.8572534603],
                                              ['C', -1.0096975537, -3.2134558984, -1.0828411727],
                                              ['C', -1.6496547071, -2.5134588595, -0.0545223130],
                                              ['C', -1.2782539691, -1.1895674728,  0.1981616803],
                                              ['C',  0.0906576267,  0.8387203882, -0.2960550684],
                                              ['C', -0.8961672869,  1.8448570513, -0.2165145863],
                                              ['C', -0.5578452911,  3.1622179253,  0.0198641969],
                                              ['C',  0.8207399895,  3.5414615134,  0.1947462265],
                                              ['C',  1.8084480415,  2.5166961566,  0.1181181691],
                                              ['C',  1.4364478511,  1.2003830241, -0.1258364018],
                                              ['O',  1.0676792904,  4.7807656709,  0.4180224177],
                                              ['O', -1.2587662420, -4.5118138191, -1.4541769692],
                                              ['C', -2.3235465711, -5.2305613788, -0.7524733671],
                                              ['C', -2.3510202348, -6.6072129004, -1.4177275690],
                                              ['C', -3.5506076805, -7.3984246184, -0.9330683249],
                                              ['S', -3.5301144993, -9.0799676721, -1.6264065105],
                                              ['O', -2.2422208883, -9.6722233426, -1.1947825968],
                                              ['O', -3.6513610781, -8.9126245785, -3.0929123728],
                                              ['O', -4.7164373593, -9.7187446365, -1.0019775156],
                                              ['H', -3.5672026693, -7.4780232701,  0.1637571944],
                                              ['H', -4.4955292697, -6.9148489369, -1.2211805879],
                                              ['H', -1.4077722137, -7.1537368289, -1.2032271649],
                                              ['H', -2.3770927024, -6.4922147170, -2.5224726281],
                                              ['H', -3.2621044660, -4.6755871706, -0.8946544955],
                                              ['H', -2.0586858440, -5.2813561955,  0.3130491969],
                                              ['H',  0.4722100115, -3.1759973440, -2.6449631741],
                                              ['H',  1.1145578311, -0.8004436405, -2.1812042020],
                                              ['H', -1.7615795694, -0.6395818795,  1.0002252855],
                                              ['H', -2.4171495115, -2.9782931891,  0.5438074496],
                                              ['H',  2.2026040693,  0.4285681014, -0.1795669520],
                                              ['H',  2.8425172867,  2.7882301550,  0.2505042335],
                                              ['H', -1.3047738539,  3.9423564809,  0.0747048789],
                                              ['H', -1.9405614027,  1.5705286241, -0.3549171054],
                                              ['C',  4.3218987075,  6.1160755576,  1.7429319942],
                                              ['C',  2.8163441920,  6.0112471044,  1.8733458824],
                                              ['C',  2.3148797604,  6.8046885644,  3.0650255974],
                                              ['C',  2.1875068722,  6.3621863379,  0.5980957539],
                                              ['Br', 2.9578193976,  8.8870044359,  0.3895573066],
                                              ['H',  2.5700158898,  6.0215360744, -0.3141517001],
                                              ['H',  1.2685227966,  6.8452510931,  0.5922196460],
                                              ['H',  4.8322661074,  5.7747030576,  2.6330102930],
                                              ['H',  4.6779214720,  5.5296111079,  0.9060544636],
                                              ['H',  4.6411069233,  7.1503593164,  1.5749950955],
                                              ['H',  2.5153272831,  4.9342114328,  2.0151022065],
                                              ['H',  2.7681390950,  6.4451204162,  3.9814771131],
                                              ['H',  2.5859314781,  7.8577559139,  2.9628064633],
                                              ['H',  1.2328828954,  6.7309549981,  3.1418868680]]
        for expected, actual in zip(expected_ts_guess_plus_1_structure, calculation_2.ts_guess_structure_plus_1):
            for index in range(1, 4):
                # TODO: Check for at least two places after the decimal point (the difference is probably due to
                # differences in the L-BFGS implementation)
                self.assertAlmostEqual(actual[index], expected[index], places=1)


if __name__ == '__main__':
    unittest.main()
