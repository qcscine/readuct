import io
import os
import subprocess
import unittest
import yaml


class Calculation:
    def __init__(self, input):
        self._input = input
        self._output = ''
        self.energy = 0
        self.lowest_frequency = ''
        self.optimized_structure = []

    def run(self):
        self.__prepare_input()
        self.__execute()
        self.__parse_output()

    def __prepare_input(self):
        with io.open('input.yaml', 'w', encoding='utf8') as outfile:
            yaml.dump(self._input, outfile, default_flow_style=False, allow_unicode=True)

    def __execute(self):
        # If the return code is non-zero, this automatically raises an exception
        self._output = subprocess.check_output(['readuct', 'input.yaml'])

    def __parse_output(self):
        splitted_output = self._output.splitlines()
        for index, line in enumerate(splitted_output):
            if line.startswith('  The (electronic) energy is'):
                self.energy = float(splitted_output[index].split()[4])
            if line == '     #    cm^-1':
               self.lowest_frequency = float(splitted_output[index + 1].split()[1])
            if self._input['tasks'][0]['type'] in ['afir', 'geoopt', 'ts']:
                self.optimized_structure = self.__load_xyz_file('default_system_output/default_system_output.xyz')

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


class TestSparrow(unittest.TestCase):
    def test_single_point_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6'
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'type': 'sp',
                            }
                          ]
                        }

        default_calculation = Calculation(default_input)
        default_calculation.run()
        self.assertAlmostEqual(default_calculation.energy, -22.314858496)

        input = { 'systems': [
                    { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                      'name': 'default_system',
                      'method': 'AM1',
                      'settings': {
                        'molecular_charge': 1,
                        'spin_multiplicity': 2
                      }
                    }
                  ],
                  'tasks': [
                    { 'input': ['default_system'],
                      'type': 'sp',
                    }
                  ]
                }

        calculation = Calculation(input)
        calculation.run()
        self.assertAlmostEqual(calculation.energy, -24.04238019)

    def test_hessian_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6',
                              'settings': {
                                'self_consistence_criterion': 1e-6
                              }
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'type': 'hessian',
                            }
                          ]
                        }

        default_calculation = Calculation(default_input)
        default_calculation.run()
        self.assertAlmostEqual(default_calculation.lowest_frequency, -156.8, places=0)

        input = { 'systems': [
                    { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                      'name': 'default_system',
                      'method': 'DFTB2',
                      'settings': {
                        'molecular_charge': -1,
                        'spin_multiplicity': 2
                      }
                    }
                  ],
                  'tasks': [
                    { 'input': ['default_system'],
                      'type': 'hessian',
                    }
                  ]
                }

        calculation = Calculation(input)
        calculation.run()
        self.assertAlmostEqual(calculation.lowest_frequency, -8039.5, places=0)

    def test_structure_optimization_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6'
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'output': ['default_system_output'],
                              'type': 'geoopt',
                            }
                          ]
                        }
        default_calculation = Calculation(default_input)
        default_calculation.run()
        expected_structure = [['O', -0.6975262492, 0.0974803942, 0.1173242033], ['O', 0.6974458952, -0.095239949, -0.1195113965], ['H', -1.0427520011, -0.5637268744, -0.5552272337], ['H', 1.0428323552, 0.5614864292, 0.5574144269]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index])

    def test_transition_state_optimization_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6',
                              'settings': {
                                'self_consistence_criterion': 1e-6
                              }
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'output': ['default_system_output'],
                              'type': 'ts',
                            },
                            { 'input': ['default_system_output'],
                              'type': 'hessian'
                            }
                          ]
                        }

        default_calculation = Calculation(default_input)
        default_calculation.run()
        expected_structure = [['O', -0.6646413836, 0.1429147063, -0.5292819464], ['O', 0.6469726001, 0.5515148483, -0.1398927703], ['H', -0.8703478537, -0.6211280378, 0.073550329], ['H', 0.8880166372, -0.0733015168, 0.5956243877]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index])
        self.assertAlmostEqual(default_calculation.lowest_frequency, -447.9, places=0)

    def test_afir_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'br_2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6',
                              'settings': {
                                'molecular_charge': -2
                              }
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'output': ['default_system_output'],
                              'type': 'afir',
                            }
                          ]
                        }

        default_calculation = Calculation(default_input)
        default_calculation.run()
        expected_structure = [['Br', -6.1090129646, 0.0, 0.0], ['Br', 6.1090129646, 0.0, 0.0]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index])

    def test_irc_task(self):
        default_input = { 'systems': [
                            { 'path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'h2o2.xyz'),
                              'name': 'default_system',
                              'method': 'PM6',
                              'settings': {
                                'self_consistence_criterion': 1e-8
                              }
                            }
                          ],
                          'tasks': [
                            { 'input': ['default_system'],
                              'output': ['default_system_ts'],
                              'type': 'ts',
                            },
                            { 'input': ['default_system_ts'],
                              'output': ['default_system_output', 'default_system_output_2'],
                              'type': 'irc',
                            },
                          ]
                        }

        default_calculation = Calculation(default_input)
        default_calculation.run()
        expected_structure = [['O', -0.410337295, 0.4861140855, -0.973816632], ['O', 0.3728130841, 1.0028853057, -0.4608172341], ['H', 0.8895162512, -0.3716114976, 1.0806016985], ['H', -0.8519920402, -1.1173878936, 0.3540321676]]
        for expected, actual in zip(expected_structure, default_calculation.optimized_structure):
            for index in range(1, 4):
                self.assertAlmostEqual(actual[index], expected[index], places=4)


if __name__ == '__main__':
    unittest.main()
