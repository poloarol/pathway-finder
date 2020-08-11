from flask import Flask, jsonify, request
from flask_cors import CORS

from flask_mail import Mail

import ast

from finder import Finder

UPLOAD_FOLDER = './static/uploads'
ALLOWED_EXT = set(['gb','gbk', 'gbff'])

# instatiation
app = Flask(__name__)
app.config.from_object(__name__)

mail = Mail(app)


# enable CORS
cors = CORS(app, resources={r'/*': {'origin': '*'}})

#sanity check
@app.route('/ping', methods=['GET'])
def ping_pong():
    return jsonify('pong!')


@app.route('/pathway', methods=['GET', 'POST'])
def index():

    """
    """

    if request.method == 'POST':
        data = list(request.form.keys())[0]
        data = ast.literal_eval(data)

        email: str = data['email'] if 'email' in data else ''
        accession: str = data['accession'] if 'accession' in data else ''
        protein: str = data['protein'] if 'protein' in data else ''
        similarity: int = data['similarity'] if 'similarity' in data else ''
        basepairs: int = data['basepairs'] if 'basepairs' in data else ''

        finder = Finder(email, accession=accession, coreGene=protein, similarity=int(similarity), bp=int(basepairs))
        paths = finder.finder()


    return paths


def allowed_file(filename: str):
    """Define the allowd files to be uploaded."""
    return filename.lower.endswith(('gb','gbk', 'gbff'))



if __name__ == '__main__':
    print("Flask App is running on port 5000")
    app.run(debug=True)
