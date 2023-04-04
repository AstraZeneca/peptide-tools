class FileExtension:
    name = 'Adam'

    def __repr__(self):
        return repr('Hello ' + self.name )

print(repr(FileExtension))