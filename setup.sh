DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PYTHONPATH="$DIR/pyext/src:$PYTHONPATH"
export PYTHONPATH

exec ${precommand} "$@"
