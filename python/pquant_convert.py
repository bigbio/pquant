#!/usr/bin/env python3

"""
This is the main tool that give access to all commands and options provided by the sparkms

@author Chakradhar Bandla

"""
import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """This is the main tool that give access to all commands and options provided by the sparkms"""


cli.add_command(json_to_parquet)
cli.add_command(psm_table)


def main():
    cli()


if __name__ == "__main__":
    main()





