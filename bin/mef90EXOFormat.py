#!/usr/bin/env python3
import click
import exodus3 as exo

@click.command()
@click.option('--preset', '-p', type=click.Choice(['vDef2D','vDef3D','None'],case_sensitive=False),default='None',help = 'preset for damage, displacement, and stress for vDef')

@click.option('--nodalvariables', '-n', type=click.STRING, help="comma separated list of variable names",default='')
@click.option('--zonalvariables', '-z', type=click.STRING, help="comma separated list of variable names",default='')
# @click.option('--filename', '-f', type=click.Path(exists=True))
@click.argument('filename', type=click.Path(exists=True))

def exoFormat(filename,nodalvariables,zonalvariables,preset):
    """This script is used to format an exodus file"""
    if preset == 'vDef2D':
        if not nodalvariables == '':
            click.echo(f"Warning: Nodal variables will be ignored for preset {preset}")
        nodalvariables = ['Damage','Displacement_X','Displacement_Y']
        if not zonalvariables == '':
            click.echo(f"Warning: Zonal variables will be ignored for preset {preset}")
        zonalvariables = ['Stress_XX','Stress_YY','Stress_XY']
    elif preset == 'vDef3D':
        if not nodalvariables == '':
            click.echo(f"Warning: Nodal variables will be ignored for preset {preset}")
        nodalvariables = ['Damage','Displacement_X','Displacement_Y','Displacement_Z']
        if not zonalvariables == '':
            click.echo(f"Warning: Zonal variables will be ignored for preset {preset}")
        zonalvariables = ['Stress_XX','Stress_YY','Stress_ZZ','Stress_YZ','Stress_XZ','Stress_XY']
    else:
        nodalvariables = nodalvariables.split(',')
        if nodalvariables == ['']:
            nodalvariables = []
        zonalvariables = zonalvariables.split(',')
        if zonalvariables == ['']:
            zonalvariables = []
    # click.echo(f"Filename: {click.format_filename(filename)}")
    # click.echo(f"Preset: {preset}")
    # click.echo(f"Nodal Variables: {nodalvariables}")
    # click.echo(f"Zonal Variables: {zonalvariables}")

    e = exo.exodus(filename, mode='a')
    e.set_global_variable_number(0)
    e.set_node_variable_number(len(nodalvariables))
    for i,name in enumerate(nodalvariables):
        e.put_node_variable_name(name,i+1)

    e.set_element_variable_number(len(zonalvariables))
    for i,name in enumerate(zonalvariables):
        e.put_element_variable_name(name,i+1)
    if len(zonalvariables) > 0:
        e.set_element_variable_truth_table([True] * e.numElemBlk.value * len(zonalvariables))

    e.close()
if __name__ == '__main__':
    exoFormat()