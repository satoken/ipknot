"""Read existing build commands and preserve ILP benchmark decoder selection."""
from functools import lru_cache
from pathlib import Path
import shlex
import subprocess


def link_command(build):
    build = Path(build)
    link = build / 'CMakeFiles/ipknot.dir/link.txt'
    if link.exists():
        return shlex.split(link.read_text())
    if (build / 'build.ninja').exists():
        commands = subprocess.check_output(
            ['ninja', '-C', str(build), '-t', 'commands', 'ipknot'], text=True)
        command = commands.splitlines()[-1]
        # CMake's Ninja link rule wraps the compiler in two no-op commands.
        command = command.removeprefix(': && ').removesuffix(' && :')
        tokens = shlex.split(command)
        if '&&' in tokens or ';' in tokens:
            raise ValueError('Unsupported compound Ninja link command')
        return tokens
    raise FileNotFoundError(f'No supported link command in {build}')


def legacy_ip_api(build):
    symbols = subprocess.check_output(
        ['nm', '-C', '-u', str(Path(build) / 'CMakeFiles/ipknot.dir/src/ipknot.cpp.o')],
        text=True)
    return 'IP::IP(IP::DirType, int, bool)' not in symbols


@lru_cache(maxsize=None)
def ilp_decoder_args(binary):
    result = subprocess.run([str(binary), '--help'], capture_output=True,
                            text=True, check=True)
    return ['--decoder', 'ilp'] if '--decoder' in result.stdout else []
