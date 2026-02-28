# Load plotly
import plotly.graph_objects as go

# Sample data
x = [1.5, 2.9, 3, 4.2, 5.6]
y = [2.2, 13.3, 4.4, 55.3, 52.1]

# Initialize a figure
fig = go.Figure()

# Add the scatter trace
fig.add_trace(go.Scatter(
    x=x, # Variable in the x-axis
    y=y, # Variable in the y-axis
    mode='markers', # This explicitly states that we want our observations to be represented by points
))

# Show
fig.show()
