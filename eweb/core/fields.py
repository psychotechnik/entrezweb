from Bio.Seq import Seq, MutableSeq

from django.db import models
from django.db.models.fields.files import FileDescriptor


class MutableSequenceField(models.TextField):

    #__metaclass__ = models.SubfieldBase

    def __init__(self, *args, **kwargs):
        super(MutableSequenceField, self).__init__(*args, **kwargs)

    def to_python(self, value):
        if isinstance(value, MutableSeq):
            return value
        return MutableSeq(value)

    def get_db_prep_value(self, value, connection, prepared=False):
        if isinstance(value, MutableSeq):
            return str(value)
        return ''

    def value_to_string(self, obj):
        value = self._get_val_from_obj(obj)
        return self.get_db_prep_value(value)

class SequenceField(models.TextField):

    #__metaclass__ = models.SubfieldBase

    def __init__(self, *args, **kwargs):
        super(SequenceField, self).__init__(*args, **kwargs)

    def to_python(self, value):
        if isinstance(value, Seq):
            return value
        return Seq(value)

    def get_db_prep_value(self, value, connection, prepared=False):
        if isinstance(value, Seq):
            return str(value)
        return ''

    def value_to_string(self, obj):
        value = self._get_val_from_obj(obj)
        return self.get_db_prep_value(value)

#from south.modelsinspector import add_introspection_rules
#add_introspection_rules([
#    ([SequenceField], [],{},),
#], ["^core\.fields\.SequenceField"])

#add_introspection_rules([
#    ([SequenceField], [],{},),
#], ["^core\.fields\.MutableSequenceField"])


#class FastaFileDescriptor(FileDescriptor):

class FastaFileField(models.FileField):
    pass

#add_introspection_rules([
#    ([SequenceField], [],{},),
#], ["^core\.fields\.FastaFileField"])



class CSVField(models.TextField):

    #__metaclass__ = models.SubfieldBase

    def __init__(self, *args, **kwargs):
        self.delimeter = kwargs.pop("delimeter", ',')
        super(CSVField, self).__init__(*args, **kwargs)

    def to_python(self, value):
        if value == '':
            return set()
        if isinstance(value, set):
            return value
        return set(map(lambda x: int(x), value.split(self.delimeter)))

    def get_db_prep_value(self, value, connection, prepared=False):
        return self.delimeter.join([str(s) for s in value])

    def value_to_string(self, obj):
        value = self._get_val_from_obj(obj)
        return self.get_db_prep_value(value)


class CSVCharField(CSVField):
    def to_python(self, value):
        if value == '':
            return set()
        if isinstance(value, set):
            return value
        return set(map(lambda x: x, value.split(self.delimeter)))

#from south.modelsinspector import add_introspection_rules
#add_introspection_rules([
#    ([CSVField], [], {},),
#], ["^core\.fields\.CSVField"])
#add_introspection_rules([
#    ([CSVField], [], {},),
#], ["^core\.fields\.CSVCharField"])

