from flask import Flask, jsonify
from flask_cors import CORS


# instatiation
app = Flask(__name__)
app.config.from_object(__name__)

# enable CORS
cors = CORS(app, resources={r'/*': {'origin': '*'}})

#sanity check
@app.route('/ping', methods=['GET'])
def ping_pong():
    return jsonify('pong!')


@app.route('/', methods=['GET'])
def index():
    return 'Hello World!!'



if __name__ == '__main__':
    print("Flask App is running on port 5000")
    app.run(debug=True)
