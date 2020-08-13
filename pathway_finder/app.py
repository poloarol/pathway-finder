from flask import Flask, jsonify, request
from flask_cors import CORS

from flask_mail import Mail

import ast
import sqlite3

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

    paths = None

    print('Request - 202 - Analysis beginning')

    if request.method == 'POST':
        data = list(request.form.keys())[0]
        data = ast.literal_eval(data)

        email: str = data['email'] if 'email' in data else ''
        accession: str = data['accession'] if 'accession' in data else ''
        protein: str = data['protein'] if 'protein' in data else ''
        similarity: int = data['similarity'] if 'similarity' in data else ''
        basepairs: int = data['basepairs'] if 'basepairs' in data else ''

        print('Launching finder ... ')

        finder = Finder(email, accession=accession, coreGene=protein, similarity=float(similarity), bp=int(basepairs))
        paths = finder.finder()

        print('Finder completed ... ')

        if len(path):
            insert_job_number(data['uuid'])
            insert_job_data(data['uuid'], paths)


        return jsonify('AMI Motivated')


@app.route('/submission', methods=['GET'])
def diagram(jobnumber):

    """
    """

    pathway = retrieve_info(job_number)
    return pathway


def allowed_file(filename: str):
    """Define the allowd files to be uploaded."""
    return filename.lower.endswith(('gb','gbk', 'gbff'))


def establish_connnection():
    
    """
    """

    try:
        sqlite3.register_adapter(uuid.UUID, lambda u: u.bytes_le)
        BASE_DIR: str = os.path.dirname(os.path.abspath(__file__))
        BASE_DIR = BASE_DIR + "/db/"
        db_path = os.path.join(BASE_DIR, "GenomeMap.db")
        connection = sqlite3.connect(db_path, check_same_thread=False)
        return connection
    except Exception:
        print(traceback.format_exc())


def insert_job_number(job_number):
    """ Insert into jobs table """
    try:
        sql = " INSERT INTO JOBS (JOB_NUMBER) VALUES (?) "
        connection = establish_connection()
        cursor = connection.cursor()
        cursor.execute(sql, (str(job_number), ))
        cursor.close()
        connection.commit()
        connection.close()
    except Exception:
        print(traceback.format_exc())

def insert_job_data(job_number, data):
    try:
        sql = " INSERT INTO GENOMEMAP_INFO (JOB_NUMBER, DATUM ) VALUES (?, ?) "
        connection = establish_connection()
        cursor = connection.cursor()
        data = json.dumps(data)
        cursor.execute(sql, (str(job_number), data))
        cursor.close()
        connection.commit()
        connection.close()
    except Exception:
        print(traceback.format_exc())

def retrieve_info(job_number):
    try:
        sql = " SELECT DATUM FROM GENOMEMAP_INFO WHERE JOB_NUMBER = ? "
        connection = establish_connection()
        cursor = connection.cursor()
        cursor.execute(sql, (str(job_number),))
        data = cursor.fetchone()
        cursor.close()
        connection.close()
        return json.loads(data[0])
    except Exception:
        print(traceback.format_exc())


if __name__ == '__main__':
    print("Flask App is running on port 5000")
    app.run(debug=True, threaded=True)
