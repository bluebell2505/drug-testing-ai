
import React from 'react';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { Button } from '@/components/ui/button';
import { Link, useNavigate } from 'react-router-dom';
import { cn } from '@/lib/utils';
import { Beaker, ArrowRight, Folder } from 'lucide-react';
import { useIsMobile } from '@/hooks/use-mobile';

const Dashboard: React.FC = () => {
  const navigate = useNavigate();
  const userName = "Dr. Smith"; // This would come from authentication in a real app
  const isMobile = useIsMobile();
  
  // Mock projects data
  const recentProjects = [
    { id: 1, name: "Project Alpha", lastUpdated: "2025-03-23", status: "In Progress" },
    { id: 2, name: "Project Beta", lastUpdated: "2025-03-22", status: "Completed" },
    { id: 3, name: "Project Gamma", lastUpdated: "2025-03-20", status: "In Progress" }
  ];

  // Mock AI stats
  const aiStats = {
    name: "DrugGen v1.0",
    status: "Ready",
    accuracy: 92,
    moleculesToday: 5,
    lastPrediction: "2025-03-24 10:30 AM"
  };

  const handleGenerateNewMolecule = () => {
    navigate('/molecule-generator');
  };

  const handleViewAllProjects = () => {
    navigate('/projects');
  };

  const handleViewAnalytics = () => {
    navigate('/analytics');
  };

  const handleViewProjectDetails = (projectId: number) => {
    navigate(`/projects/${projectId}`);
  };

  const handleGenerateMoleculeForProject = (projectId: number) => {
    navigate(`/molecule-generator?projectId=${projectId}`);
  };
  
  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6 overflow-y-auto">
          {/* Welcome Section */}
          <div className="bg-gray-800/50 border border-gray-700 rounded-xl p-6 mb-6">
            <h1 className="text-white text-3xl font-bold mb-4">Welcome back, {userName}!</h1>
            <div className="flex flex-wrap gap-6 mb-6">
              <div className="text-blue-400">3 Active Projects</div>
              <div className="text-blue-400">5 Molecules Generated Today</div>
              <div className="text-blue-400">AI Model: Ready</div>
            </div>
            <Button 
              className="w-full sm:w-auto bg-blue-600 hover:bg-blue-700 text-white flex items-center justify-center gap-2 p-3 h-auto rounded-lg text-base font-medium"
              onClick={handleGenerateNewMolecule}
            >
              <Folder className="w-5 h-5" />
              Generate New Molecule
            </Button>
          </div>

          {/* Recent Projects Section */}
          <h2 className="text-white text-2xl font-bold mb-4">Recent Projects</h2>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
            {recentProjects.map(project => (
              <div 
                key={project.id} 
                className="bg-gray-800/50 border border-gray-700 p-5 rounded-xl"
              >
                <h3 className="text-white font-medium text-xl mb-2">{project.name}</h3>
                <p className="text-gray-400 text-sm mb-1">Last updated: {project.lastUpdated}</p>
                <p className={cn(
                  "text-sm font-medium mb-4",
                  project.status === "In Progress" ? "text-blue-400" : 
                  project.status === "Completed" ? "text-green-400" : "text-yellow-400"
                )}>
                  Status: {project.status}
                </p>
                <div className="flex flex-col sm:flex-row gap-2">
                  <Button 
                    variant="outline" 
                    className="text-blue-400 border-blue-400 hover:bg-blue-900/20 text-sm h-9" 
                    onClick={() => handleViewProjectDetails(project.id)}
                  >
                    View Details
                  </Button>
                  <Button 
                    className="bg-blue-600 hover:bg-blue-700 text-white text-sm h-9" 
                    onClick={() => handleGenerateMoleculeForProject(project.id)}
                  >
                    Generate Molecule
                  </Button>
                </div>
              </div>
            ))}
          </div>
          <Button 
            variant="ghost" 
            className="text-blue-400 hover:bg-blue-900/20 mb-6 flex items-center gap-1"
            onClick={handleViewAllProjects}
          >
            View All Projects <ArrowRight size={14} />
          </Button>

          {/* AI Model Insights Section */}
          <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl mb-6">
            <h2 className="text-white text-2xl font-bold mb-4">AI Model Insights</h2>
            <div className="flex flex-col">
              <div className="flex justify-between mb-2">
                <span className="text-white">Model Name:</span>
                <span className="text-blue-400">{aiStats.name}</span>
              </div>
              <div className="flex justify-between mb-2">
                <span className="text-white">Status:</span>
                <span className="text-blue-400">{aiStats.status}</span>
              </div>
              <div className="flex justify-between mb-2">
                <span className="text-white">Prediction Accuracy:</span>
                <span className="text-blue-400">{aiStats.accuracy}%</span>
              </div>
              <div className="flex justify-between mb-2">
                <span className="text-white">Molecules Generated Today:</span>
                <span className="text-blue-400">{aiStats.moleculesToday}</span>
              </div>
              <div className="flex justify-between mb-4">
                <span className="text-white">Last Prediction:</span>
                <span className="text-blue-400">{aiStats.lastPrediction}</span>
              </div>
            </div>
            <Button 
              variant="ghost" 
              className="text-blue-400 hover:bg-blue-900/20 flex items-center gap-1"
              onClick={handleViewAnalytics}
            >
              View Analytics <ArrowRight size={14} />
            </Button>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default Dashboard;
