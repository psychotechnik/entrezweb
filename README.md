# Flask API Project Management Application

This project is a **reference architecture** for building Flask applications. It demonstrates a clean, maintainable structure for a Flask app where users can create projects and tasks.


## Table of Contents

- [Flask API Project Management Application](#flask-api-project-management-application)
  - [Table of Contents](#table-of-contents)
  - [Project Overview](#project-overview)
  - [Testing Strategy](#testing-strategy)
    - [Unit Tests](#unit-tests)
    - [Integration Tests](#integration-tests)
  - [Application Structure](#application-structure)
    - [Main File](#main-file)
    - [Router Files](#router-files)
    - [Service Layer](#service-layer)
    - [Domain Logic](#domain-logic)
  - [Installation and running](#installation-and-running)
  - [Running Tests](#running-tests)

## Project Overview

The application allows users to create projects add tasks.

This project serves as a **reference architecture** to demonstrate best practices for organizing Flask projects, focusing on separation of concerns, testing strategies, and modularization. 

## Testing Strategy

The goal is to achieve **high test coverage** at every layer, ensuring that each component works in isolation.

****
### Unit Tests

- Each layer of the architecture (services, routers, etc.) is tested individually.
- We only test **public methods** to avoid coupling tests to internal implementation details, 
- No external service calls are made during unit tests.

### Integration Tests
- Integration testing will be handled by the front end's end-to-end testing.
- In the case where we can not perform front end tests, use a **few integration tests**  to test the entire flow with real API calls and database interactions.
- These tests are minimal and only serve as sanity checks to ensure the app functions as expected in a real-world scenario.
- Integration tests ensure that the different layers of the architecture work together correctly.

## Application Structure

The application is designed to be modular, with clear separation of concerns. Here's a breakdown of the main components:

### Main File

The `projects.py` file serves as the entry point for the application. It initializes the Flask app, makes any app-wide configurations, and includes the necessary routers:

```python
from projects.entrypoints.flask import create_app 
app = create_app()
```

### Router Files

The routers are responsible for defining the API routes. They remain very clean and only handles request validation and forwarding the call to the service layer:

```python
from projects.adapters.users import repository
from projects.entrypoints.flask import db
from projects.entrypoints.flask.api.auth import token_auth
from projects.service_layer.users import handlers

@bp.route('/users/<int:id>', methods=['GET'])
@token_auth.login_required
def get_user(id):
    repo = repository.FlaskSqlAlchemyRepository(db)
    user = handlers.get_user(id, repo)
    return user.to_dict()
```

### Service Layer

The service layer acts as the **main coordinator** of business logic. It handles calls to the connected repositories such as the FlaskSqlAlchemyRepository:

```python
def get_user_by_username(username: str, repo: repository.AbstractRepository):
    return repo.get_by_username(username)
```

### Domain Logic

The domain logic is responsible for handling business rules. The domain logic ensures that the core business rules are applied consistently across the application.

For example, the `user.py` file handles user domain model:

```python
class User:
    id: Optional[int]
    password_hash: Optional[str]
    token: Optional[str]
    token_expiration: Optional[datetime]

    def __init__(self, 
            username: str, 
            email: str, 
            is_manager: bool = False,
                ) -> None:

        self.username = username
        self.email = email
        self.is_manager = is_manager

    def __repr__(self):
        return '<User {}>'.format(self.username)

```

## Installation and running

To run the project locally, follow these steps:

1. Clone the repository:

2. Create and activate a virtual environment:

   - On macOS/Linux:

     ```bash
     python -m venv .venv
     source .venv/bin/activate
     ```

   - On Windows:

     ```bash
     python -m venv .venv
     .venv\Scripts\activate
     ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Set up the environment variables:

   - Copy the `example.flaskenv` file and rename it to `.flaskenv`.
   - Edit the `.flaskenv` file to include your PostgreSQL connection information:

   ```
   DATABASE_URL=postgresql://<username>@localhost/projecs_app
   ```

5. Set up the database (using Alembic with flask-migrate for migrations):

   ```bash
   flask db upgrade
   ```

6. Start the Flask application:

   ```bash
   flask --app src/projects/entrypoints/flask/projects run --debug
   ```

You can now visit `http://127.0.0.1:8034/swagger` to interact with the API through the automatically generated Swagger UI.

## Running Tests

To run all tests, use the following command:

```bash
pytest
```

To run only the **unit tests**:

```bash
pytest tests/unit_tests/
```

To run the **integration tests**:

```bash
pytest tests/integration_tests/
```
