#!/usr/bin/env python3
"""Optional certificate generation using Singular; verification does not need it."""
import gzip,json,shutil,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    executable=shutil.which('Singular')
    if executable is None:raise SystemExit('Install Singular (or use the supplied verified certificate).')
    data=json.loads((ROOT/'data/f210-equations.json').read_text())
    output=ROOT/'generated-certificates';output.mkdir(exist_ok=True)
    script='ring r=0,('+','.join(data['variables'])+'),dp;\n'
    script+='ideal E='+','.join(data['reduced_equations'])+';\n'
    script+='option(redSB);matrix T;ideal G=liftstd(E,T);\n'
    script+='write("basis.txt",string(G));write("lift.txt",string(T));quit;\n'
    with tempfile.TemporaryDirectory(prefix='tpe-certificate-') as temporary:
        tmp=Path(temporary);(tmp/'generate.sing').write_text(script)
        result=subprocess.run([executable,'-q','generate.sing'],cwd=tmp,text=True,capture_output=True,check=True)
        if '?' in result.stdout or '?' in result.stderr:raise RuntimeError(result.stdout+result.stderr)
        shutil.copyfile(tmp/'basis.txt',output/'f210-basis.txt')
        with (tmp/'lift.txt').open('rb') as source,(output/'f210-lift.txt.gz').open('wb') as target:
            with gzip.GzipFile(filename='',mode='wb',fileobj=target,mtime=0,compresslevel=9) as compressed:
                shutil.copyfileobj(source,compressed)
    print('Generated certificates:',output)
    print('They are separate from the frozen certificates and SHA256 manifest.')

if __name__=='__main__':main()
