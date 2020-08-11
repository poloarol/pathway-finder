""" Flask forms """

from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField, SubmitField, RadioField
from wtforms.validators import DataRequired
from flask_wtf.file import FileField, FileRequired, FileAllowed

class InfoForm(FlaskForm):
    """
    """

    email = StringField('E-mail', validators=[DataRequired('Provide your e-mail address')])
    accession_number = StringField('Accession Number')
    gene_of_interest = StringField('Provide the sequence of interest')
    similarity = IntegerField("% identity tolerated", default=0.65)
    basepairs = IntegerField("Number of bases", default=2500)
    upload = FileField('File Upload', validators=[FileAllowed(['gb', 'gbk', 'gbff'], 'GenBank Only')])
    submit = SubmitField('Submit')