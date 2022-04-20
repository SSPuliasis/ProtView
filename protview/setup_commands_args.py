import argparse
import protview.setup_commands_main as setup_commands_main

parser = argparse.ArgumentParser(
    description="adds user-defined proteases to the RPG protease options", add_help=False)

def main():
    setup_commands_main.add_ud_proteases()
    print('user-defined proteases added')