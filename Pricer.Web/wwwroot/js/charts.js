// wwwroot/js/charts.js
// Render the binomial-tree convergence chart via Chart.js

let convergenceChart = null;

window.renderConvergenceChart = function (canvasId, labels, binomData, bsData) {
    const ctx = document.getElementById(canvasId);
    if (!ctx) return;

    if (convergenceChart) {
        convergenceChart.data.labels = labels;
        convergenceChart.data.datasets[0].data = binomData;
        convergenceChart.data.datasets[1].data = bsData;
        convergenceChart.update();
        return;
    }

    convergenceChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [
                {
                    label: 'Binomial tree',
                    data: binomData,
                    borderColor: '#0d6efd',
                    backgroundColor: 'transparent',
                    pointRadius: 4,
                    tension: 0.25,
                    borderWidth: 2,
                },
                {
                    label: 'Black-Scholes',
                    data: bsData,
                    borderColor: '#198754',
                    backgroundColor: 'transparent',
                    borderDash: [8, 4],
                    pointRadius: 0,
                    borderWidth: 2,
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { position: 'top' },
                tooltip: {
                    callbacks: {
                        label: ctx => `${ctx.dataset.label}: $${ctx.parsed.y.toFixed(4)}`
                    }
                }
            },
            scales: {
                x: { title: { display: true, text: 'Number of steps (n)' } },
                y: {
                    title: { display: true, text: 'Option price ($)' },
                    ticks: { callback: v => '$' + v.toFixed(2) }
                }
            }
        }
    });
};
