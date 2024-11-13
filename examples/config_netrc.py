import os
import getpass

netrcpath = os.path.expanduser('~/.netrc')
cookiespath = os.path.expanduser('~/.urs_cookies')
dodsrcpath = os.path.expanduser('~/.dodsrc')

# Make empty files
if os.path.exists(netrcpath):
    print(netrcpath + ' exists; delete to remake')
else:
    with open(netrcpath, 'w') as outf:
        outf.write('')
    os.chmod(netrcpath, 0x600)
    with open(netrcpath, 'w') as outf:
        outf.write(f"""
machine urs.earthdata.nasa.gov
  login {getpass.getpass('URS User:')}
  password {getpass.getpass('URS Password:')}
""")

if os.path.exists(cookiespath):
    print(cookiespath + ' exists; delete to remake')
else:
    with open(cookiespath, 'w') as outf:
        outf.write('')

if os.path.exists(dodsrcpath):
    print(dodsrcpath + ' exists; delete to remake')
else:
    with open(dodsrcpath, 'w') as outf:
        outf.write(f"""HTTP.NETRC={netrcpath}
HTTP.COOKIEJAR={cookiespath}
""")
