#!/usr/bin/env python3
"""Compare overlapping animal-host and sugar-bait interventions with an old binary.

Usage: python3 test/feeding_regression.py --build build --original PATH --output NEW_DIR
The output directory must be new. Uses the existing larviciding scenario as input.
"""
import argparse
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--original', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=False)
source = Path(__file__).resolve().parent
tree = ET.parse(source / 'scenarioNoMPDLarviciding.xml')
root = tree.getroot()
interventions = root.find('interventions')
for species in interventions.findall('vectorPop/intervention/description/anopheles'):
    for field, initial in [('probAdditionalDeathSugarFeeding', '0.12'),
                           ('seekingDeathRateIncrease', '0.23'), ('probDeathOvipositing', '0.17')]:
        effect = ET.SubElement(species, field, initial=initial)
        ET.SubElement(effect, 'decay', L='0.5', function='exponential', CV='0.2')
changes = ET.SubElement(interventions, 'nonHumanHostsModifications')
for index, initial in enumerate(['0.13', '0.27']):
    intervention = ET.SubElement(changes, 'intervention', name=f'animal-effect-{index}',
                                nonHumanHostsName='unprotectedAnimals')
    ET.SubElement(intervention, 'decay', L='0.5', function='exponential', CV='0.4')
    description = ET.SubElement(intervention, 'description')
    for species in root.findall('entomology/vector/anopheles'):
        animal = ET.SubElement(description, 'anopheles', mosquito=species.attrib['mosquito'])
        for field in ['availabilityReduction', 'preprandialKillingEffect',
                      'postprandialKillingEffect', 'restingKillingEffect', 'fecundityReduction']:
            ET.SubElement(animal, field, initial=initial)
    timed = ET.SubElement(intervention, 'timed')
    for date in ['0', '155', '465']:
        ET.SubElement(timed, 'deploy', time=date)
scenario = args.output.resolve() / 'feeding.xml'
tree.write(scenario, encoding='UTF-8', xml_declaration=True)
for output_format in ['bin', 'text']:
    baseline = args.output / ('original-' + output_format)
    command = [sys.executable, str(source / 'exact_regression.py'), '--build', str(args.build),
               '--format', output_format, str(scenario)]
    subprocess.run(command + ['--executable', str(args.original), '--output', str(baseline)], check=True)
    subprocess.run(command + ['--baseline', str(baseline), '--output',
                   str(args.output / ('current-' + output_format))], check=True)
