class PyMetabolismError(StandardError):
    def __init__(self, msg="", *args, **kwargs):
        super(PyMetabolismError, self).__init__(*args, **kwargs)
#        self.errno
        self.strerror = msg

    def __str__(self):
        print self.strerror