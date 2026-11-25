# -*- coding: utf-8 -*-
"""Update TDCRPy after its window has closed, then reopen it.

Run by :mod:`tdcrpy.gui` from a COPY in a temporary folder, never from the
package, and never imported: the point is to replace the package, so nothing
here may import ``tdcrpy`` or live inside it while pip works.

Why the window cannot simply run pip itself
-------------------------------------------
On Windows a file that a process holds open cannot be removed. The window's
launcher (``tdcrpy-gui.exe``) is such a file while the window is open, and
``TDCR_model_lib`` keeps its decay-data archives open for as long as any
process -- a calculation, a notebook kernel, a script -- has it imported.
pip then fails halfway through the uninstall with ``WinError 32`` and cannot
roll back cleanly: it leaves a ``~dcrpy`` folder and a mixed installation
behind (observed with 2.20.35 -> 2.20.36). So this helper

1. waits for the window's process to end, then for its launcher to be
   released;
2. checks that no OTHER program still holds the package's files, and stops
   before pip touches anything if one does;
3. runs pip, checks the version it installed, says what happened, and
   reopens the window.
"""
import argparse
import datetime
import glob
import os
import subprocess
import sys
import time

_NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)


def wait_for_pid(pid, timeout=120.0):
    """Wait until process *pid* has exited. Returns ``False`` on timeout."""
    if os.name == 'nt':
        import ctypes
        SYNCHRONIZE = 0x00100000
        k32 = ctypes.windll.kernel32
        handle = k32.OpenProcess(SYNCHRONIZE, False, int(pid))
        if not handle:
            return True                          # already gone
        try:
            return k32.WaitForSingleObject(handle, int(timeout * 1000)) == 0
        finally:
            k32.CloseHandle(handle)
    end = time.time() + timeout
    while time.time() < end:
        try:
            os.kill(int(pid), 0)
        except OSError:
            return True
        time.sleep(0.2)
    return False


def is_locked(path):
    """Whether another process holds *path* so that pip could not remove
    it. Tested by renaming it away and back: a rename is exactly what pip's
    uninstall does first, and it fails on the same lock."""
    if os.name != 'nt' or not os.path.exists(path):
        return False
    probe = path + '.tdcrpy-update-check'
    try:
        os.rename(path, probe)
    except OSError:
        return True
    os.rename(probe, path)
    return False


def wait_until_free(paths, timeout=30.0):
    """Wait for *paths* to be released; return those still locked."""
    end = time.time() + timeout
    while True:
        busy = [p for p in paths if is_locked(p)]
        if not busy or time.time() > end:
            return busy
        time.sleep(0.5)


def package_files(package_dir):
    """The files of the installed package most likely to be held open."""
    if not package_dir:
        return []
    return sorted(glob.glob(os.path.join(package_dir, 'decayData', '*.zip')))


def installed_version(python):
    try:
        out = subprocess.run(
            [python, '-c', 'from importlib.metadata import version; '
                           'print(version("TDCRPy"))'],
            capture_output=True, text=True, timeout=120,
            creationflags=_NO_WINDOW)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def tell(title, text, error=False, dialog=True):
    """Say how it went: the window is closed, so a dialog is the only
    channel the user will see."""
    if not dialog:
        print(f'{title}: {text}')
        return
    try:
        import tkinter as tk
        from tkinter import messagebox
        root = tk.Tk()
        root.withdraw()
        (messagebox.showerror if error else messagebox.showinfo)(title, text)
        root.destroy()
    except Exception:                                 # no Tk: nothing to do
        print(f'{title}: {text}')


def relaunch(python):
    """Reopen the window without a console behind it."""
    gui_python = python
    if os.name == 'nt':
        w = os.path.join(os.path.dirname(python), 'pythonw.exe')
        if os.path.exists(w):
            gui_python = w
    flags = 0
    if os.name == 'nt':
        flags = (getattr(subprocess, 'DETACHED_PROCESS', 0)
                 | getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0))
    subprocess.Popen([gui_python, '-m', 'tdcrpy.gui'], close_fds=True,
                     creationflags=flags, cwd=os.path.expanduser('~'))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--pid', type=int, required=True,
                    help='process of the window to wait for')
    ap.add_argument('--python', required=True,
                    help='interpreter of the TDCRPy installation')
    ap.add_argument('--spec', default='TDCRPy',
                    help='what to install, e.g. TDCRPy==2.20.37')
    ap.add_argument('--from-version', default='?')
    ap.add_argument('--package-dir', default='')
    ap.add_argument('--launcher', default='',
                    help='the tdcrpy-gui launcher, released after the window')
    ap.add_argument('--log', required=True)
    ap.add_argument('--relaunch', action='store_true')
    ap.add_argument('--no-dialog', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args(argv)
    dialog = not args.no_dialog
    title = 'TDCRPy update'

    def log(text):
        with open(args.log, 'a', encoding='utf-8') as fh:
            fh.write(f'[{datetime.datetime.now():%H:%M:%S}] {text}\n')

    log(f'updating {args.from_version} -> {args.spec} with {args.python}')
    if not wait_for_pid(args.pid):
        log('the window did not close; giving up')
        tell(title, 'The TDCRPy window did not close, so the update was not '
                    'started. Nothing was changed.', True, dialog)
        return 1

    busy = wait_until_free([args.launcher] if args.launcher else [])
    busy += wait_until_free(package_files(args.package_dir), timeout=5.0)
    if busy:
        log('still in use: ' + ', '.join(busy))
        tell(title,
             'TDCRPy is still in use by another program, so the update was '
             'not started and nothing was changed.\n\nClose any Jupyter '
             'notebook, Python console or script that has imported tdcrpy '
             '(Windows locks the files it reads), then use Help ▸ Check for '
             'updates again.\n\nIn use:\n' + '\n'.join(busy), True, dialog)
        if args.relaunch:
            relaunch(args.python)
        return 2

    cmd = [args.python, '-m', 'pip', 'install', '--upgrade',
           '--disable-pip-version-check', '--no-input', args.spec]
    log('running: ' + ' '.join(cmd))
    with open(args.log, 'a', encoding='utf-8') as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT,
                             creationflags=_NO_WINDOW)
    now = installed_version(args.python)
    log(f'pip exited with {rc}; installed version now {now}')

    if rc == 0 and now:
        tell(title, f'TDCRPy {args.from_version} was updated to {now}.'
                    + ('\n\nThe window will now reopen.' if args.relaunch
                       else ''), False, dialog)
        if args.relaunch:
            relaunch(args.python)
        return 0
    tell(title, f'The update did not complete (pip exit code {rc}; installed '
                f'version now {now or "unknown"}).\n\nThe full pip output is '
                f'in:\n{args.log}\n\nTo finish by hand, close every program '
                f'using TDCRPy and run:\n'
                f'"{args.python}" -m pip install --upgrade {args.spec}',
         True, dialog)
    if args.relaunch and now:
        relaunch(args.python)
    return 1


if __name__ == '__main__':
    sys.exit(main())
