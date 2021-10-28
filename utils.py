import os
import subprocess


class ApplicationError(subprocess.CalledProcessError):
    """Raised when an application returns a non-zero exit status (OBSOLETE).
    The exit status will be stored in the returncode attribute, similarly
    the command line string used in the cmd attribute, and (if captured)
    stdout and stderr as strings.
    This exception is a subclass of subprocess.CalledProcessError.
    >>> err = ApplicationError(-11, "helloworld", "", "Some error text")
    >>> err.returncode, err.cmd, err.stdout, err.stderr
    (-11, 'helloworld', '', 'Some error text')
    >>> print(err)
    Non-zero return code -11 from 'helloworld', message 'Some error text'
    """

    def __init__(self, returncode, cmd, stdout="", stderr=""):
        """Initialize the class."""
        self.returncode = returncode
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr

    def __str__(self):
        """Format the error as a string."""
        # get first line of any stderr message
        try:
            msg = self.stderr.lstrip().split("\n", 1)[0].rstrip()
        except Exception:  # TODO, ValueError? AttributeError?
            msg = ""
        if msg:
            return "Non-zero return code %d from %r, message %r" % (
                self.returncode,
                self.cmd,
                msg,
            )
        else:
            return "Non-zero return code %d from %r" % (self.returncode, self.cmd)

    def __repr__(self):
        """Represent the error as a string."""
        return "ApplicationError(%i, %s, %s, %s)" % (
            self.returncode,
            self.cmd,
            self.stdout,
            self.stderr,
        )


def run_cmd(cmd, stdin=None, stdout=True, stderr=True, cwd=None, env=None):
    """Execute command, wait for it to finish, return (stdout, stderr).
    Runs the command line tool and waits for it to finish. If it returns
    a non-zero error level, an exception is raised. Otherwise two strings
    are returned containing stdout and stderr.
    The optional stdin argument should be a string of data which will be
    passed to the tool as standard input.
    The optional stdout and stderr argument may be filenames (string),
    but otherwise are treated as a booleans, and control if the output
    should be captured as strings (True, default), or ignored by sending
    it to /dev/null to avoid wasting memory (False). If sent to a file
    or ignored, then empty string(s) are returned.
    The optional cwd argument is a string giving the working directory
    to run the command from. See Python's subprocess module documentation
    for more details.
    The optional env argument is a dictionary setting the environment
    variables to be used in the new process. By default the current
    process' environment variables are used. See Python's subprocess
    module documentation for more details.
    Default example usage::
        from Bio.Emboss.Applications import WaterCommandline
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
                                     stdout=True, auto=True,
                                     asequence="a.fasta", bsequence="b.fasta")
        print("About to run: %s" % water_cmd)
        std_output, err_output = water_cmd()
    This functionality is similar to subprocess.check_output(). In general
    if you require more control over running the command, use subprocess
    directly.
    When the program called returns a non-zero error level, a custom
    ApplicationError exception is raised. This includes any stdout and
    stderr strings captured as attributes of the exception object, since
    they may be useful for diagnosing what went wrong.
    """
    if not stdout:
        stdout_arg = open(os.devnull, "w")
    elif isinstance(stdout, str):
        stdout_arg = open(stdout, "w")
    else:
        stdout_arg = subprocess.PIPE

    if not stderr:
        stderr_arg = open(os.devnull, "w")
    elif isinstance(stderr, str):
        if stdout == stderr:
            stderr_arg = stdout_arg  # Write both to the same file
        else:
            stderr_arg = open(stderr, "w")
    else:
        stderr_arg = subprocess.PIPE

    child_process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=stdout_arg,
        stderr=stderr_arg,
        universal_newlines=True,
        cwd=cwd,
        env=env,
        shell=True,
    )
    # Use .communicate as can get deadlocks with .wait(), see Bug 2804
    stdout_str, stderr_str = child_process.communicate(stdin)
    if not stdout:
        assert not stdout_str, stdout_str
    if not stderr:
        assert not stderr_str, stderr_str
    return_code = child_process.returncode

    # Particularly important to close handles on Jython and PyPy
    # (where garbage collection is less predictable) and on Windows
    # (where cannot delete files with an open handle):
    if not stdout or isinstance(stdout, str):
        # We opened /dev/null or a file
        stdout_arg.close()
    if not stderr or (isinstance(stderr, str) and stdout != stderr):
        # We opened /dev/null or a file
        stderr_arg.close()

    if return_code:
        raise ApplicationError(return_code, cmd, stdout_str, stderr_str)
    return stdout_str, stderr_str
