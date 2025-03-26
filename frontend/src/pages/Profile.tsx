
import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { 
  Avatar, 
  AvatarFallback, 
  AvatarImage 
} from '@/components/ui/avatar';
import { 
  Button,
  Card,
  CardHeader,
  CardContent,
  CardTitle,
  CardDescription,
  CardFooter,
  Alert,
  AlertTitle,
  AlertDescription
} from '@/components/ui';
import { 
  BarChart3, 
  Medal, 
  Clock, 
  Edit, 
  Share2, 
  PieChart, 
  UserPlus,
  Calendar,
  Lightbulb,
  Settings,
  Award,
  Bookmark,
  CheckCircle,
  Target
} from 'lucide-react';
import { useToast } from '@/hooks/use-toast';
import { EditProfileModal } from '@/components/EditProfileModal';

// Mock data for the profile
const profileData = {
  user: {
    name: "Dr. John Smith",
    email: "john.smith@example.com",
    institution: "Research Institute of Molecular Science",
    role: "Lead Researcher",
    joinedDate: "March 1, 2025",
    initials: "JS"
  },
  projects: [
    { id: 1, name: "Project Alpha", stage: "In Testing", molecules: 12, lastUpdated: "2025-03-23" },
    { id: 2, name: "Project Beta", stage: "Completed", molecules: 8, lastUpdated: "2025-03-22" },
    { id: 3, name: "Project Gamma", stage: "Preclinical Validation", molecules: 5, lastUpdated: "2025-03-20" }
  ],
  achievements: [
    { id: 1, name: "Molecule Master", description: "Generated 50+ molecules", icon: "Medal", completed: true },
    { id: 2, name: "Project Pioneer", description: "Contributed to 5+ projects", icon: "Medal", completed: true },
    { id: 3, name: "Binding Expert", description: "Achieved high binding affinity for 10 molecules", icon: "Medal", completed: false },
    { id: 4, name: "Toxicity Tamer", description: "Reduced toxicity to Low for 5 molecules", icon: "Medal", completed: false }
  ],
  recentActivity: [
    { id: 1, description: "Generated a molecule for Project Beta", date: "2025-03-22", link: "/projects/2" },
    { id: 2, description: "Optimized a molecule in Project Alpha", date: "2025-03-20", link: "/projects/1" },
    { id: 3, description: "Completed Project Gamma", date: "2025-03-15", link: "/projects/3" }
  ],
  insights: [
    { id: 1, title: "Most Active Target Protein", content: "HER2 (45% of molecules target HER2)" },
    { id: 2, title: "Success Rate", content: "85% of generated molecules have high binding affinity" },
    { id: 3, title: "Recommendation", content: "Consider optimizing molecules for EGFR to diversify your portfolio" }
  ],
  preferences: {
    targetProteins: ["EGFR", "HER2"],
    optimizationGoals: ["High Solubility", "Low Toxicity"],
    notificationSettings: "Email notifications for project updates: Enabled"
  },
  progressMetrics: {
    contribution: 75, // percentage
    projectDistribution: [
      { stage: "In Testing", percentage: 40 },
      { stage: "Completed", percentage: 30 },
      { stage: "On Hold", percentage: 20 },
      { stage: "Planning", percentage: 10 }
    ]
  },
  collaboration: {
    collaborators: [
      { id: 1, name: "Dr. Jane Doe", project: "Project Alpha" },
      { id: 2, name: "Dr. Michael Chen", project: "Project Beta" }
    ]
  },
  learningResources: [
    { id: 1, title: "Learn more about HER2 binding affinity optimization", link: "#" },
    { id: 2, title: "Advanced techniques for reducing toxicity", link: "#" }
  ],
  gamification: {
    rank: 3,
    challenge: "Generate 5 molecules with high solubility this week to earn the Solubility Star badge"
  },
  visualization: {
    favoriteMolecule: {
      smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
      bindingAffinity: "High for EGFR (Score: 9.2)",
      name: "Aspirin derivative"
    }
  },
  customization: {
    bio: "Passionate about AI-driven drug discovery and targeting HER2 for cancer therapies"
  }
};

const Profile: React.FC = () => {
  const navigate = useNavigate();
  const { toast } = useToast();
  const [editModalOpen, setEditModalOpen] = useState(false);

  // Share profile function
  const handleShareProfile = () => {
    const profileUrl = `${window.location.origin}/profile/123`;
    navigator.clipboard.writeText(profileUrl);
    toast({
      title: "Success",
      description: "Profile link copied to clipboard."
    });
  };

  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6 overflow-y-auto">
          {/* Header */}
          <div className="flex justify-between items-center mb-8">
            <div className="flex items-center gap-4">
              <Avatar className="h-16 w-16 bg-blue-500 text-white">
                <AvatarFallback className="text-xl font-bold">{profileData.user.initials}</AvatarFallback>
              </Avatar>
              <div>
                <h1 className="text-white text-3xl font-bold">User Profile</h1>
                <p className="text-gray-400">View and manage your profile information</p>
              </div>
            </div>
            <div className="flex gap-3">
              <Button 
                variant="outline" 
                className="text-white border-gray-700"
                onClick={handleShareProfile}
              >
                <Share2 size={18} className="mr-2" />
                Share Profile
              </Button>
              <Button 
                className="bg-blue-500 hover:bg-blue-600 text-white"
                onClick={() => setEditModalOpen(true)}
              >
                <Edit size={18} className="mr-2" />
                Edit Profile
              </Button>
            </div>
          </div>
          
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* Left Column */}
            <div className="space-y-6 lg:col-span-1">
              {/* Personal Information */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader>
                  <CardTitle className="text-white">Personal Information</CardTitle>
                </CardHeader>
                <CardContent className="text-gray-300 space-y-3">
                  <div>
                    <p className="text-gray-400 text-sm">Full Name</p>
                    <p className="font-medium">{profileData.user.name}</p>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Email</p>
                    <p>{profileData.user.email}</p>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Institution</p>
                    <p>{profileData.user.institution}</p>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Role</p>
                    <p>{profileData.user.role}</p>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Joined Date</p>
                    <p>{profileData.user.joinedDate}</p>
                  </div>
                </CardContent>
              </Card>
              
              {/* Achievements */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <CardTitle className="text-white">Achievements</CardTitle>
                </CardHeader>
                <CardContent className="pt-2">
                  <div className="space-y-3">
                    {profileData.achievements.map(achievement => (
                      <div 
                        key={achievement.id} 
                        className={`p-3 rounded-lg flex items-center gap-3 ${
                          achievement.completed ? 'bg-blue-900/30' : 'bg-gray-700/30'
                        }`}
                      >
                        <div className={`p-2 rounded-full ${
                          achievement.completed ? 'bg-blue-500' : 'bg-gray-600'
                        }`}>
                          <Medal size={18} className="text-white" />
                        </div>
                        <div className="flex-1">
                          <p className="text-white font-medium">{achievement.name}</p>
                          <p className="text-gray-400 text-sm">{achievement.description}</p>
                        </div>
                        {achievement.completed && (
                          <CheckCircle size={18} className="text-blue-400" />
                        )}
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
              
              {/* Preferences */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <div className="flex justify-between items-center">
                    <CardTitle className="text-white">Preferences</CardTitle>
                    <Button variant="ghost" size="icon" className="text-gray-400 hover:text-white">
                      <Settings size={18} />
                    </Button>
                  </div>
                </CardHeader>
                <CardContent className="pt-2 text-gray-300 space-y-3">
                  <div>
                    <p className="text-gray-400 text-sm">Preferred Target Proteins</p>
                    <div className="flex flex-wrap gap-2 mt-1">
                      {profileData.preferences.targetProteins.map((protein, idx) => (
                        <span key={idx} className="px-2 py-1 bg-blue-900/30 rounded text-sm text-blue-300">
                          {protein}
                        </span>
                      ))}
                    </div>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Default Optimization Goals</p>
                    <div className="flex flex-wrap gap-2 mt-1">
                      {profileData.preferences.optimizationGoals.map((goal, idx) => (
                        <span key={idx} className="px-2 py-1 bg-green-900/30 rounded text-sm text-green-300">
                          {goal}
                        </span>
                      ))}
                    </div>
                  </div>
                  <div>
                    <p className="text-gray-400 text-sm">Notification Settings</p>
                    <p className="text-sm">{profileData.preferences.notificationSettings}</p>
                  </div>
                </CardContent>
              </Card>
              
              {/* Bio */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <CardTitle className="text-white">Bio</CardTitle>
                </CardHeader>
                <CardContent className="pt-2 text-gray-300">
                  <p className="text-sm italic">{profileData.customization.bio}</p>
                </CardContent>
              </Card>
            </div>
            
            {/* Center Column */}
            <div className="space-y-6 lg:col-span-2">
              {/* Progress Metrics */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <CardTitle className="text-white">Progress Metrics</CardTitle>
                </CardHeader>
                <CardContent className="pt-2">
                  <div className="mb-4">
                    <p className="text-gray-400 text-sm mb-2">Drug Discovery Progress</p>
                    <div className="h-4 bg-gray-700 rounded-full overflow-hidden">
                      <div 
                        className="h-full bg-blue-500 rounded-full" 
                        style={{ width: `${profileData.progressMetrics.contribution}%` }}
                      ></div>
                    </div>
                    <p className="text-gray-400 text-sm mt-1">{profileData.progressMetrics.contribution}% towards Preclinical Candidate Selection</p>
                  </div>
                  
                  <div>
                    <p className="text-gray-400 text-sm mb-3">Project Distribution</p>
                    <div className="grid grid-cols-2 gap-2">
                      {profileData.progressMetrics.projectDistribution.map((item, idx) => (
                        <div key={idx} className="flex items-center gap-2">
                          <div className={`w-3 h-3 rounded-full ${
                            idx === 0 ? 'bg-blue-500' : 
                            idx === 1 ? 'bg-green-500' : 
                            idx === 2 ? 'bg-yellow-500' : 'bg-purple-500'
                          }`}></div>
                          <span className="text-gray-300 text-sm">{item.stage}: {item.percentage}%</span>
                        </div>
                      ))}
                    </div>
                  </div>
                </CardContent>
              </Card>
              
              {/* Projects Overview */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <CardTitle className="text-white">Projects Worked On</CardTitle>
                </CardHeader>
                <CardContent className="pt-2">
                  <div className="space-y-4">
                    {profileData.projects.map(project => (
                      <div key={project.id} className="bg-gray-700/50 rounded-lg p-4">
                        <div className="flex justify-between mb-2">
                          <h3 className="text-white font-medium">{project.name}</h3>
                          <span className={`text-sm px-2 py-0.5 rounded-full ${
                            project.stage === "In Testing" ? "bg-blue-500/20 text-blue-300" :
                            project.stage === "Completed" ? "bg-green-500/20 text-green-300" :
                            project.stage === "Preclinical Validation" ? "bg-purple-500/20 text-purple-300" : 
                            "bg-yellow-500/20 text-yellow-300"
                          }`}>
                            {project.stage}
                          </span>
                        </div>
                        <div className="flex justify-between text-gray-400 text-sm mb-3">
                          <span>Molecules: {project.molecules}</span>
                          <span>Last Updated: {project.lastUpdated}</span>
                        </div>
                        <Button 
                          variant="outline" 
                          className="w-full text-blue-400 border-blue-500/50 hover:bg-blue-900/20"
                          onClick={() => navigate(`/projects/${project.id}`)}
                        >
                          View Project
                        </Button>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
              
              {/* Activity Timeline and AI Insights */}
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Activity Timeline */}
                <Card className="bg-gray-800 border-gray-700">
                  <CardHeader className="pb-2">
                    <CardTitle className="text-white">Recent Activity</CardTitle>
                  </CardHeader>
                  <CardContent className="pt-2">
                    <div className="space-y-4 relative before:absolute before:left-3 before:top-0 before:h-full before:w-[1px] before:bg-gray-700">
                      {profileData.recentActivity.map(activity => (
                        <div key={activity.id} className="pl-8 relative">
                          <div className="absolute left-0 top-1 h-6 w-6 rounded-full flex items-center justify-center bg-blue-500/20 border border-blue-500">
                            <Clock size={12} className="text-blue-500" />
                          </div>
                          <p className="text-white">{activity.description}</p>
                          <p className="text-gray-400 text-sm">{activity.date}</p>
                          <Button 
                            variant="link" 
                            className="text-blue-400 px-0 hover:text-blue-300"
                            onClick={() => navigate(activity.link)}
                          >
                            View Details
                          </Button>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
                
                {/* AI Insights */}
                <Card className="bg-gray-800 border-gray-700">
                  <CardHeader className="pb-2">
                    <CardTitle className="text-white">AI-Powered Insights</CardTitle>
                  </CardHeader>
                  <CardContent className="pt-2">
                    <div className="space-y-4">
                      {profileData.insights.map(insight => (
                        <div key={insight.id} className="p-3 bg-blue-900/10 rounded-lg border border-blue-500/20">
                          <p className="text-blue-400 font-medium mb-1">{insight.title}</p>
                          <p className="text-gray-300 text-sm">{insight.content}</p>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
              </div>
              
              {/* Collaboration and Learning Resources */}
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Collaborators */}
                <Card className="bg-gray-800 border-gray-700">
                  <CardHeader className="pb-2">
                    <CardTitle className="text-white">Collaborators</CardTitle>
                  </CardHeader>
                  <CardContent className="pt-2">
                    <div className="space-y-3">
                      {profileData.collaboration.collaborators.map(collaborator => (
                        <div key={collaborator.id} className="flex items-center justify-between p-3 bg-gray-700/30 rounded-lg">
                          <div className="flex items-center gap-3">
                            <Avatar className="h-8 w-8 bg-gray-600 text-white">
                              <AvatarFallback>{collaborator.name.charAt(0)}</AvatarFallback>
                            </Avatar>
                            <div>
                              <p className="text-white">{collaborator.name}</p>
                              <p className="text-gray-400 text-xs">Worked on {collaborator.project}</p>
                            </div>
                          </div>
                          <Button variant="ghost" size="icon" className="text-gray-400 hover:text-white">
                            <UserPlus size={18} />
                          </Button>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
                
                {/* Learning Resources */}
                <Card className="bg-gray-800 border-gray-700">
                  <CardHeader className="pb-2">
                    <CardTitle className="text-white">Recommended Resources</CardTitle>
                  </CardHeader>
                  <CardContent className="pt-2">
                    <div className="space-y-3">
                      {profileData.learningResources.map(resource => (
                        <a 
                          key={resource.id} 
                          href={resource.link}
                          className="flex items-center gap-3 p-3 bg-gray-700/30 rounded-lg hover:bg-gray-700/50 transition"
                        >
                          <div className="p-2 bg-blue-500/20 rounded-full">
                            <Lightbulb size={16} className="text-blue-400" />
                          </div>
                          <p className="text-blue-300 text-sm">{resource.title}</p>
                        </a>
                      ))}
                    </div>
                  </CardContent>
                </Card>
              </div>
              
              {/* Gamification */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <div className="flex justify-between items-center">
                    <CardTitle className="text-white">Achievements & Challenges</CardTitle>
                    <p className="text-amber-400 font-medium">Rank #{profileData.gamification.rank}</p>
                  </div>
                </CardHeader>
                <CardContent className="pt-2">
                  <Alert className="bg-amber-900/20 border-amber-500/30 text-amber-300">
                    <Award className="h-4 w-4" />
                    <AlertTitle>Current Challenge</AlertTitle>
                    <AlertDescription>
                      {profileData.gamification.challenge}
                    </AlertDescription>
                  </Alert>
                </CardContent>
              </Card>
              
              {/* Favorite Molecule */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader className="pb-2">
                  <CardTitle className="text-white">Favorite Molecule</CardTitle>
                </CardHeader>
                <CardContent className="pt-2">
                  <div className="flex flex-col md:flex-row gap-4">
                    <div className="bg-gray-700 rounded-md p-4 h-40 flex items-center justify-center flex-1">
                      <p className="text-gray-400 text-center">Molecule visualization placeholder</p>
                    </div>
                    <div className="flex-1">
                      <h3 className="text-white font-medium mb-2">{profileData.visualization.favoriteMolecule.name}</h3>
                      <p className="text-gray-400 text-sm mb-1">SMILES</p>
                      <p className="text-gray-300 text-sm break-all mb-3">{profileData.visualization.favoriteMolecule.smiles}</p>
                      <p className="text-gray-400 text-sm mb-1">Binding Affinity</p>
                      <p className="text-gray-300">{profileData.visualization.favoriteMolecule.bindingAffinity}</p>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </div>
          </div>
        </main>
      </div>
      
      {/* Edit Profile Modal */}
      <EditProfileModal 
        isOpen={editModalOpen}
        onClose={() => setEditModalOpen(false)}
        profileData={profileData}
      />
      
      <Footer />
    </div>
  );
};

export default Profile;
