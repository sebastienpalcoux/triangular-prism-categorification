#!/usr/bin/env python3
"""Apply spectrum criteria or generate localization equations for new rings."""
import argparse
import json
from pathlib import Path
import sys
from fusion_ring import load_fusion_ring, validate_fusion_ring
from spectrum import scan


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Exact criteria and localization equations for an arbitrary fusion ring.',
        epilog='Input: N[i][j][k], unit index 0; raw JSON tensor or object with fusion_matrices. '
               'An object may include labels and positive rational dimensions. Noncommutative rings are accepted. '
               'A search without witnesses does not prove categorifiability.')
    commands = parser.add_subparsers(dest='command', required=True)
    spectrum = commands.add_parser('spectrum', help='search the zero/one-spectrum criteria')
    spectrum.add_argument('input', type=Path, help='fusion-rule JSON file')
    spectrum.add_argument('--criterion', choices=('zero','one','both'), default='both')
    spectrum.add_argument('--mode', choices=('first','count','all'), default='count',
                        help='first: stop after each requested obstruction is found; '
                             'count: exhaustive counts and first witnesses (default); '
                             'all: also retain all witnesses')
    spectrum.add_argument('--output', type=Path, help='write the JSON report here instead of stdout')
    localize = commands.add_parser('localize', help='generate exact necessary localization equations')
    localize.add_argument('input', type=Path, help='fusion-rule JSON file, optionally with localization settings')
    localize.add_argument('--config', type=Path, help='separate localization-specification JSON; overrides input settings')
    localize.add_argument('--output', type=Path, help='write the polynomial-system JSON here instead of stdout')
    localize.add_argument('--singular', type=Path, help='also export the system as a Singular input file')
    args = parser.parse_args(argv)
    try:
        N, metadata = load_fusion_ring(args.input)
        validation = validate_fusion_ring(N, metadata.get('dimensions'))
        if args.command == 'spectrum':
            result = scan(N, args.criterion, args.mode, validate=False)
        else:
            from localization import generate_from_spec, to_singular
            spec = json.loads(args.config.read_text()) if args.config else metadata.get('localization')
            if not isinstance(spec, dict):
                raise ValueError('Provide a localization object in the input or --config SPEC.json')
            result = generate_from_spec(N, spec)
            if args.singular:
                if args.output and args.output.resolve() == args.singular.resolve():
                    raise ValueError('--output and --singular must name different files')
                args.singular.write_text(to_singular(result))
        result['validation'] = validation
        for key in ('name','labels'):
            if key in metadata:
                result[key] = metadata[key]
        output = json.dumps(result, indent=2)+'\n'
        if args.output:
            args.output.write_text(output)
        else:
            sys.stdout.write(output)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
